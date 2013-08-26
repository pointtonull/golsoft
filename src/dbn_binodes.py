"""
This file contains everything needed for the binet version of the DBN,
based on Pietros non-binet DBN node implementation.
"""

import numpy as np

import mdp
import bimdp
from dbn_nodes import DBNLayerNode


class DBNLayerBiNode(bimdp.BiNode, DBNLayerNode):
    """Adapter to turn the DBNLayerNode into a BiNode."""

    def __init__(self, node_id, hidden_dim, visible_dim=None, dtype=None):
        super(DBNLayerBiNode, self).__init__(node_id=node_id,
                                             hidden_dim=hidden_dim,
                                             visible_dim=visible_dim,
                                             dtype = dtype)
    
    def _is_bi_learning(self):
        return self._updown_initialized

    # v as the first argument is filled with the x value
    def _up_pass(self, v, epsilon=0.1, decay=0.0, momentum=0.0):
        v, pv, deltav = super(DBNLayerBiNode, self)._up_pass(
                                                v, epsilon, decay, momentum)
        # pv and deltav are currently not used,
        # but could be passed in the message
        return v
    
    # h has a different dimension than x, so have to transport it in message,
    # x is None, but is always given as the first argument
    def _down_pass(self, x, h, top_updates=0, epsilon=0.1, decay=0.0,
                   momentum=0.0):
        v, pv, deltav = super(DBNLayerBiNode, self)._down_pass(
                                    h, top_updates, epsilon, decay, momentum)
        # pv and deltav are currently not used
        return None, {"h": h}
        

class DBNMasterBiNode(bimdp.BiNode):
    """Node sits atop the DBN and manages the updown training phase."""
    
    def __init__(self, dbn_ids, sender_id, node_id="dbn_master",
                 input_dim=None, output_dim=None, dtype=None):
        """Initialize node.
        
        sender_id -- id of the sender node at the flow bottom.
        dbn_ids -- List with the ids of all the DBN node in the correct order.
        """
        self.dbn_ids = dbn_ids
        self.sender_id = sender_id
        super(DBNMasterBiNode, self).__init__(node_id=node_id,
                                              input_dim=input_dim,
                                              output_dim=input_dim,
                                              dtype=dtype)
    
    @bimdp.binode_coroutine(["msg_x", "msg"])
    def _train(self, x, msg_x, max_iter, min_error, msg):
        """Manage the DBN training."""
        i = 0
        error = np.inf
        while i < max_iter and error > min_error:
            ## up execution phase
            orig_x = msg_x
            msg = {}
            for dbn_id in self.dbn_ids:
                msg[dbn_id + "->method"] = "up_pass"
            x, _, _ = yield orig_x, msg, self.sender_id
            ## down execution phase
            msg = {}
            for dbn_id in self.dbn_ids:
                msg[dbn_id + "->target"] = -1
                msg[dbn_id + "->method"] = "down_pass"
            msg[self.sender_id + "->target"] = self.node_id
            msg[self.sender_id + "->no_x"] = True  # avoid x dimension error
            msg["h"] = x
            x, _, _ = yield None, msg, -1
            ## execution phase
            # do one normal execution for the error calculation
            x, msg_x, _ = yield orig_x, None, self.sender_id
            ## inverse phase
            msg = {}
            for dbn_id in self.dbn_ids:
                msg[dbn_id + "->method"] = "inverse"
            msg[self.sender_id + "->target"] = self.node_id
            msg[self.sender_id + "->no_x"] = True
            x, _, _ = yield x, msg, -1
            ## calculate new error and restart up phase
            i += 1
            error = float(mdp.numx.absolute(orig_x - msg_x).sum())
        ## this should end the training
        raise StopIteration()


@mdp.extension_method("html", DBNMasterBiNode, "_html_representation")
def master_html_representation(self):
    if self._coroutine_instances and "_train" in self._coroutine_instances:
        co_locals = self._coroutine_instances["_train"].gi_frame.f_locals
        return (['iter counter: %d' % co_locals["i"],
                 'error: %.5f' % co_locals["error"]])
    else:
        return ""

def get_DBN_flow(n_layers, hidden_dims):
    """Factory function for DBNs."""
    dbn_ids = []
    nodes = [bimdp.nodes.SenderBiNode(node_id="sender")]
    for i_layer in range(n_layers):
        dbn_ids.append("dbn_%d" % (i_layer+1))
        nodes.append(DBNLayerBiNode(node_id=dbn_ids[i_layer],
                                    hidden_dim=hidden_dims[i_layer]))
    nodes.append(DBNMasterBiNode(dbn_ids=dbn_ids,
                                 sender_id="sender",
                                 node_id="dbn_master"))
    return bimdp.BiFlow(nodes)   
