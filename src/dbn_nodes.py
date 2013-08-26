"""
This file contains the original DBN node from Pietro.
"""

import mdp
from mdp import numx
from mdp.utils import mult

random = mdp.numx_rand.random
exp = mdp.numx.exp

# TODO: DBNLayer with labels
# TODO: DBNLayer with Gaussian visible variables
# TODO: remove self._rbm after greedy phase
# TODO: some functionality of RBMNode is duplicated with small variation;
#    should this be a subclass of RBMNode?

class DBNLayerNode(mdp.Node):

    def __init__(self, hidden_dim, visible_dim = None, dtype = None):
        super(DBNLayerNode, self).__init__(input_dim = visible_dim,
                                           output_dim = hidden_dim,
                                           dtype = dtype)
        # delegates computations to an RBM during the greedy phase
        self._rbm = mdp.nodes.RBMNode(hidden_dim, visible_dim, dtype)
        self._updown_initialized = False

    def sample_h(self, v):
        """Sample the hidden variables given observations v.

        Returns a tuple (prob_h, h), where prob_h[n,i] is the
        probability that variable 'i' is one given the observations
        v[n,:], and h[n,i] is a sample from the posterior probability."""
        self._pre_execution_checks(v)
        return self._sample_h(v)

    def sample_v(self, h):
        """Sample the observed variables given hidden variable state h.

        Returns a tuple (prob_v, v), where prob_v[n,i] is the
        probability that variable 'i' is one given the hidden variables
        h[n,:], and v[n,i] is a sample from that conditional probability."""
        self._pre_inversion_checks(h)
        return self._sample_v(h)

    def _sample_h(self, v):
        # P(h=1|v,W,b)
        probs = 1./(1. + exp(-self.bh - mult(v, self.w_rec)))
        h = (probs > random(probs.shape)).astype('d')
        return probs, h

    def _sample_v(self, h):
        # P(v=1|h,W,b)
        probs = 1./(1. + exp(-self.bv - mult(h, self.w_gen.T)))
        v = (probs > random(probs.shape)).astype('d')
        return probs, v
    
    def _execute(self, v, return_probs=True):
        """If 'return_probs' is True, returns the probability of the
        hidden variables h[n,i] being 1 given the observations v[n,:].
        If 'return_probs' is False, return a sample from that probability.
        """
        probs, h = self._sample_h(v)
        if return_probs:
            return probs
        else:
            return h

    def _inverse(self, h):
        probs, v = self._sample_v(h)
        return v

    # greedy phase training: delegate to RBM
    def _train(self, x, epsilon=0.1, decay=0., momentum=0.):
        self._rbm.train(x, epsilon=epsilon, decay=decay, momentum=momentum)

    def _stop_training(self):
        self._rbm.stop_training()
        self._init_updown()
    # /greedy phase training

    def _init_updown(self):
        """
        The greedy phase of learning is delegated to an RBM.
        Here the learned weights are decoupled into generative and
        recognition weights, and other parameters are initialized.
        """
        self._updown_initialized = True
        # recognition and generative weights
        self.w_rec = self._rbm.w.copy()
        self.w_gen = self._rbm.w.copy()
        # biases
        self.bv = self._rbm.bv.copy()
        self.bh = self._rbm.bh.copy()
        # changes in parameters during learning
        # used to add momentum
        self.dw_wake = 0.
        self.dw_sleep = 0.
        self.dbv = self.dbh = 0.

    # this corresponds to the wake phse
    def _up_pass(self, v, epsilon=0.1, decay=0., momentum=0.):
        """
        Returns (sample from hidden layer, norm of the weight change)
        """
        self._pre_execution_checks(v)

        # sample from hidden layer
        ph, h = self._sample_h(v)
        # reconstruct input
        pv1, v1 = self.sample_v(h)
        
        # adapt generative weights
        delta = mult((v - pv1).T, h)/v.shape[0]
        self.dw_wake = (momentum*self.dw_wake
                        + epsilon*(delta - decay*self.w_gen))
        self.w_gen += self.dw_wake

        # adapt biases
        delta = (v - pv1).mean(axis=0)
        self.dbv = momentum*self.dbv + epsilon*delta
        self.bv += self.dbv

        return h, ph, mdp.utils.norm2(self.dbv)

    # this corresponds to the sleep phase
    def _down_pass(self, h, top_updates=0, epsilon=0.1, decay=0., momentum=0.):
        """
        top_updates -- set >0 for top node, so that it ends up sampling
                       from the prior
        """
        # TODO: check input

        pv, v = self._sample_v(h)
        for _ in range(top_updates):
            ph, h = self._sample_h(v)
            pv, v = self._sample_v(h)
            
        # reconstruct hidden state
        ph1, h1 = self._sample_h(v)
        
        # adapt generative weights
        delta = mult(v.T, (h - ph1))/v.shape[0]
        self.dw_sleep = (momentum*self.dw_sleep
                         + epsilon*(delta - decay*self.w_rec))
        self.w_rec += self.dw_sleep

        # adapt biases
        delta = (h - ph1).mean(axis=0)
        self.dbh = momentum*self.dbh + epsilon*delta
        self.bh += self.dbh
        
        return v, pv, mdp.utils.norm2(self.dbh)
