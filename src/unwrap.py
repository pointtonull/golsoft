#!/usr/bin/env python
#-*- coding: UTF-8 -*-

from bisect import insort

from blist import blist
import cv2
import image
import numpy as np

pi = np.pi
tau = 2 * pi


def sdct(image):
    """
    Secure discrete cosine trasform
    """
    try:
        return cv2.dct(image)
    except:
        return cv2.dct(image[:-1, :-1])


def unwrap_wls(wrapped):
    """
    Weighted least squares algoritm.
    The fastest one but extremly innacurate.
    TODO: implement weights!!
    """
    rows, cols = wrapped.shape

    rowdiff = np.concatenate((np.diff(wrapped, 1, 0), np.zeros((1, cols))), 0)
    coldiff = np.concatenate((np.diff(wrapped, 1, 1), np.zeros((rows, 1))), 1)

    wrowdiff = np.mod(rowdiff + pi, tau) - pi
    wcoldiff = np.mod(coldiff + pi, tau) - pi

    rhox = np.diff(np.concatenate((np.zeros((1, rows)), wrowdiff), 0), axis=0)
    rhoy = np.diff(np.concatenate((np.zeros((cols, 1)), wcoldiff), 1), axis=1)

    rho = rhox + rhoy
    dct = sdct(rho)

    col = np.mgrid[pi / cols:pi + pi / cols: pi / cols]
    row = np.mgrid[pi / rows:pi + pi / rows: pi / rows]
    cosines = 2 * (np.cos(row)[:, np.newaxis] + np.cos(col) - 2)

    try:
        phiinv = dct / cosines
    except:
        phiinv = dct / cosines[:-1, :-1]
    unwrapped = cv2.idct(phiinv)

    return unwrapped


def unwrap_qg(phase, quality_map, equalize=True, margin=5):
    """
    Quality Guided Path Following unwrapping algoritm
    This algoritm uses the correlation array as quality map to guide the
    unwrapping path avoiding the tricky zones.

    Note: Correlation as also know as module image.

    Returns the unwrapped phase.
    """
    
    assert phase.shape == quality_map.shape
    shape = phase.shape
    rows, cols = shape

    phase /= tau

    if equalize:
        quality_map = image.equalize(quality_map)

    unwrappeds = np.zeros_like(phase).astype(bool)
    maxpos = np.unravel_index(quality_map.argmax(), quality_map.shape)
    unwrappeds[maxpos] = True
    while not unwrappeds.all():
        to_unwrap = binary_dilation(unwrappeds) - unwrappeds
        threshold = (quality_map * to_unwrap).max() - margin
        to_unwrap = quality_map >= threshold
        showimage(normalize(phase), normalize(unwrappeds.astype(int)))
        phase, unwrappeds = unwrap_mask(phase, unwrappeds, to_unwrap)
    return phase


"""
Implementarlo ahora se sale del cronograma pero creo que tengo una idea que
captura la misma información que multi-path quality guided unwrap y es muy
vectorizable.

    1. Crear un mapa de pendientes a partir del mapa de calidad. Este mapa de
       pendientes codifica de dirección en que se debe derivar e integrar. En cada
       pixel codifica guarda un número del [1, 4] donde 1 es arriba, 2, derecha,
       ect. La dirección a escoger debe ser siempre la que cree la de mejor
       calidad (el vecino con mayor valor en el mapa de calidad).
    2. Crear la derivada discreta siguiendo el mapa de pendientes.
    3. Modificar la derivada para que todos los números de acerquen lo más posible
       a 0 eliminando las congruensias mayores. matriz =% pi ¿?
    4. Integrar 3. siguiendo el mapa de pendientes.

Lo truculento está en la integración. Todos los pasos previos son sencillos y
relativamente faciles de programar. Al integrar de puede partir del máximo
absoluto, asignarle 0. Luego integrar sus vecinos en orden de peor a mejor
calidad (para evitar que los pixeles con mala calidad sobreescriban a los que
tienen buena calidad). Repetir hasta haber integrado todo el mapa. Las áreas del
mapa con mala calidad se aislan solas rompiendo el avance uniforme pero caen al ser
rodeadas.
"""

#def get_branch_cuts(residue_charge, max_box_radius, im_mask):

#    """
#    Generates branch cuts based on the phase residues. This is done using
#    the Goldstein method, as described in "Two-dimensional phase unwrapping:
#    theory, algorithms and software" by Dennis Ghiglia and Mark Pritt.
#    "residue_charge" is a matrix wherein positive residues are 1 and negative
#    residues are 0.  "max_box_radius" defines the maximum search radius for the
#    balancing of residues. If this is too large, areas will be isolated by the
#    branch cuts.  "IM_mask" is a binary matrix. This serves as an artificial
#    border for the branch cuts to connect to.
#    """

#    rowdim, coldim = residue_charge.shape

#    branch_cuts = np.zeros_like(IM_mask)
#    # Define initial branch cuts borders as the mask.
#    residue_charge_masked = residue_charge
#    residue_charge[IM_mask != 0] = 0
#    # Remove all residues except those in the mask
#    cluster_counter = 1
#    # Keep track of the number of residues in each cluster
#    satellite_residues = 0
#    # Keep track of the number of satellite residues accounted for

#    residue_binary = residue_charge != 0
#    # Logical matrix indicating the position of the residues
#    residue_balanced = np.zeros_like(IM_mask)
#    # Initially designate all residues as unbalanced

#    rowres, colres = residue_binary.nonzero()
#    # Find the coordinates of the residues
#    adjacent_residues = np.zeros(rowdim, coldim)
#    # Defines the positions of additional residues found in the search box
#    missed_residues = 0
#    # Keep track of the effective number of residues left unbalanced because of

#    print('Calculating branch cuts ...')
#    temp = rowres.shape
#    for i in xrange(1, temp(0) + 1):
#        # Loop through the residues
#        radius = 1
#        # Set the initial box size
#        r_active = rowres[i]
#        # Coordinates of the active residue
#        c_active = colres[i]
#        count_nearby_residues_flag = 1
#        # Flag to indicate whether or not to keep track of the nearby residues
#        cluster_counter = 1
#        # Reset the cluster counter
#        adjacent_residues = np.zeros(rowdim, coldim)
#        # Reset the adjacent residues indicator
#        charge_counter = residue_charge_masked[r_active, c_active]
#        # Store the initial residue charge
#        if residue_balanced[r_active, c_active] != 1:
#            # Has this residue already been balanced?
#            while charge_counter != 0: # Loop until balanced
#                # This portion of code serves to search the box perimeter,
#                # place branch cuts, and keep track of the summed residue charge
#                for m in range(r_active - radius, r_active +radius + 1):
#                    # Coordinates of the box border pixels
#                    # *** I COULD MAKE THIS MORE EFFICIENT***
#                    for n in range(c_active - radius, c_active + radius + 1):
#                        if (abs(m - r_active) == radius or abs(n - c_active) == radius) and charge_counter != 0:
#                            # Ensure that only the border pixels are being scrutinised
#                            if m <= 1 or m >= rowdim or n <= 1 or n >= coldim:
#                                # Is the current pixel on the image border?
#                                if m >= rowdim:
#                                    m = rowdim
#                                    # Make sure that the indices aren't too large for the matrix
#                                if n > coldim:
#                                    n = coldim
#                                if n < 1:
#                                    n = 1
#                                if m < 1:
#                                    m = 1
#                                branch_cuts = PlaceBranchCutInternal(branch_cuts, r_active, c_active, m, n)
#                                # Place a branch cut between the active point and the border
#                                cluster_counter += 1
#                                # Keep track of how many residues have been clustered
#                                charge_counter = 0
#                                # Label the charge as balanced
#                                residue_balanced(r_active, c_active) = 1
#                                # Mark the centre residue as balanced
#                            if IM_mask(m, n) == 0:
#                                branch_cuts = PlaceBranchCutInternal(branch_cuts,
#                                    r_active, c_active, m, n)
#                                    # Place a branch cut between the active point and the mask border
#                                cluster_counter += 1
#                                # Label the charge as balanced
#                                residue_balanced(r_active, c_active) = 1
#                                # Mark the centre residue as balanced
#                            if residue_binary[m, n] == 1:
#                                # Is the current pixel a residue?
#                                if count_nearby_residues_flag == 1:
#                                    # Are we keeping track of the residues encountered?
#                                    adjacent_residues[m, n] = 1
#                                    # Mark the current residue as adjacent
#                                branch_cuts = PlaceBranchCutInternal(branch_cuts,
#                                    r_active, c_active, m, n)
#                                # Place a branch cut regardless of the charge_counter value
#                                cluster_counter += 1
#                                # Keep track of how many residues have been clustered
#                                if residue_balanced[m, n] == 0:
#                                    residue_balanced[m, n] = 1
#                                    # Mark the current residue as balanced
#                                    charge_counter += residue_charge_masked[m, n]
#                                    #Update the charge counter
#                                if charge_counter == 0:
#                                    # Is the active residue balanced?
#                                    residue_balanced(r_active, c_active) = 1
#                                    # Mark the centre (active) residue as balanced

#                # This next portion of code centres the box on the adjacent
#                # residues. If the charge is still not balanced after moving
#                # through all adjacent residues, increase the box radius and
#                # centre the box around the initial residue.

#                if sum(sum(adjacent_residues)) == 0:
#                    # If there are no adjacent residues:
#                    radius += 1
#                    # Enlarge the box
#                    r_active = rowres[i]
#                    # Centre the larger box about the original active residue
#                    c_active = colres[i]
#                else:
#                    # If there are adjacent residues:
#                    if count_nearby_residues_flag == 1:
#                        # Run this bit once per box being searched
#                        r_adjacent, c_adjacent = adjacent_residues.nonzero()
#                        # Find the coordinates of the adjacent residues
#                        adjacent_size = r_adjacent.shape
#                        # %How many residues are on the perimeter
#                        r_active = r_adjacent[1]
#                        # Centre the search box about a nearby residue
#                        c_active = c_adjacent[1]
#                        adjacent_residue_count = 1
#                        residue_balanced[r_active, c_active] = 1
#                        # Mark the centre (active) residue as balanced before moving on
#                        count_nearby_residues_flag = 0
#                        # Disable further counting of adjacent residues
#                    else:
#                        adjacent_residue_count += 1
#                        # Move to the next nearby residue
#                        if adjacent_residue_count<=adjacent_size[1]:
#                            r_active = r_adjacent[adjacent_residue_count]
#                            # Centre the search box about the next nearby residue
#                            c_active = c_adjacent[adjacent_residue_count]
#                        else:
#                            # Ok, we're done with this box
#                            radius += 1
#                            # Enlarge the box and move on
#                            r_active = rowres[i]
#                            # Centre the larger box about the original active residue
#                            c_active = colres[i]
#                            adjacent_residues = np.zeros(rowdim, coldim)
#                            # Reset the adjacent residues indicator
#                            count_nearby_residues_flag = 1
##                           Enable further counting of adjacent residues
#                if radius > max_box_radius:
##                   Enforce a maximum boundary condition
#                    if cluster_counter != 1:
##                       The single satellite residues will still be taken care of
#                        missed_residues = missed_residues + 1
##                       This effectively means that we have an unbalanced residue
#                    else:
#                        satellite_residues += 1
#                        # This residues is about to be accounted for...
#                    charge_counter = 0

##                   This next portion of code ensures that single satellite
##                   residues are not left unaccounted for. The box is simply
##                   expanded regardless until the border or ANY other residue
##                   is found.

#                    while cluster_counter == 1:
##                        If the centre residue is a single satellite, then
##                        continue to increase the search radius
#                         r_active = rowres[i]
##                        Centre the box on the original residue
#                         c_active=colres[i]
#                         for m in range(r_active-radius, r_active + radius + 1):
##                           Coordinates of the box border pixels
##                           ***This code works but could definitely be made more
##                           efficient***
#                            for n in range(c_active - radius, c_active + radius + 1):
#                                if (abs(m - r_active) == radius
#                                    or abs(n - c_active) == radius):
##                                   Ensure that only the border pixels are
##                                   being scrutinised
#                                    if m <= 1 or m >= rowdim or n <= 1 or n >= coldim:
##                                       Is the current pixel on the image border?
#                                        if m > rowdim:
#                                            m = rowdim
##                                           Make sure that the indices aren't too
##                                           large for the matrix
#                                        if n > coldim:
#                                            n = coldim
#                                        if n < 1:
#                                            n = 1
#                                        if m < 1:
#                                            m = 1
#                                        branch_cuts = PlaceBranchCutInternal(
#                                            branch_cuts, r_active, c_active, m, n)
##                                       Place a branch cut between the active
##                                       point and the border
#                                        cluster_counter += 1
#                                        residue_balanced[r_active, c_active] = 1
##                                       Mark the centre (active) residue as balanced
#                                    if IM_mask(m, n) == 0:
##                                       Does the point fall on the mask
#                                        branch_cuts = PlaceBranchCutInternal(
#                                            branch_cuts, r_active, c_active, m, n)
##                                       Place a branch cut between the active
##                                       point and the mask border
#                                        cluster_counter += 1
##                                       Keep track of how many residues have
##                                       been clustered
#                                        residue_balanced[r_active, c_active] = 1
##                                       Mark the centre (active) residue as balanced
#                                    if residue_binary(m,n) == 1:
##                                       Is the current pixel a residue?
#                                        branch_cuts = PlaceBranchCutInternal(
#                                            branch_cuts, r_active, c_active, m, n)
##                                       Place a branch cut regardless of the
##                                       charge_counter value
#                                        cluster_counter += 1
##                                       Keep track of how many residues have been
##                                       clustered
#                                        residue_balanced[r_active, c_active] = 1
##                                       Mark the centre (active) residue as balanced
#                         radius += 1
##                        Enlarge the box and continue searching


#def get_branch_cuts(branch_cuts, r1, c1, r2, c2):
#    """
#    PlaceBranchCutInternal.m places a branch cut between the points [r1, c1] and
#    [r2, c2]. The matrix branch_cuts is binary, with 1's depicting a
#    branch cut.
#    """

#    branch_cuts[r1, c1] = 1
##   Fill the initial points
#    branch_cuts[r2, c2] = 1
#    radius = sqrt((r2 - r1) ** 2 + (c2 - c1) ** 2)
##   Distance between points
#    m = (c2 - c1)/(r2 - r1)
##   Line gradient
#    theta = atan(m)
##   Line angle

#    if c1 == c2:
##       Cater for the infinite gradient
#        if r2 > r1:
#            for i in range(1, radius + 1):
#                r_fill = r1 + i
#                branch_cuts[r_fill, c1] = 1
#            return branch_cuts
#        else:
#            for i in range(1, radius + 1):
#                r_fill = r2 + i
#                branch_cuts[r_fill, c1] = 1
#            return branch_cuts

##   Use different plot functions for the different quadrants (This is very clumsy,
##   I'm sure there is a better way of doing it...)

#    if theta < 0 and c2 < c1:
#        for i in range(1, radius + 1):
##           Number of points to fill in
#            r_fill = r1 + round(i * cos(theta))
#            c_fill = c1 + round(i * sin(theta))
#            branch_cuts[r_fill, c_fill] = 1
#    elif theta < 0 and c2 > c1:
#         for i in range(1, radius + 1):
##           Number of points to fill in
#            r_fill = r2 + round(i * cos(theta))
#            c_fill = c2 + round(i * sin(theta))
#            branch_cuts[r_fill, c_fill] = 1
#    elif theta > 0 and c2 > c1:
#        for i in range(1, radius + 1):
##           Number of points to fill in
#            r_fill = r1 + round(i * cos(theta))
#            c_fill = c1 + round(i * sin(theta))
#            branch_cuts(r_fill, c_fill) = 1
#    elif theta > 0 and c2 < c1:
#         for i in range(1, radius + 1):
##           Number of points to fill in
#            r_fill = r2 + round(i * cos(theta))
#            c_fill = c2 + round(i * sin(theta))
#            branch_cuts[r_fill, c_fill] = 1
#    else:
#        raise NotImplementedError("Must never happend")
#    return branch_cuts


#def FloodFill(IM_phase, branch_cuts, IM_mask):
#    """
#    returns:
#        IM_unwrapped,
#        rowref,
#        colref
#    FloodFill.m unwraps the phase image, avoiding all branch cuts.
#    """

#    r_dim, c_dim = IM_phase.shape
#    xref, yref = ginput[1] # escogido punto de referencia
#    colref = round(xref)
#    rowref = round(yref)

#    if branch_cuts(rowref,colref) == 1:
#        raise ValueError('Selected point corresponds to a branch cut.')

#    IM_unwrapped = np.zeros((r_dim, c_dim))
#    unwrapped_binary = np.zeros((r_dim, c_dim))
#    adjoin = zeros((r_dim, c_dim))

#    adjoin[rowref - 1, colref] = 1
##   Label the first four adjoining pixels
#    adjoin[rowref + 1, colref] = 1
#    adjoin[rowref, colref - 1] = 1
#    adjoin[rowref, colref + 1] = 1
#    IM_unwrapped[rowref, colref] = IM_phase[rowref, colref]
##   Mark the first pixel as unwrapped
#    unwrapped_binary[rowref, colref] = 1

#    print('Performing floodfill operation ...')
#    count_limit = 0
#    adjoin_stuck = 0
#    while sum(adjoin[range(2, r_dim), range(2, c_dim)]) != 0:
##       Loop until there are no adjoining pixels or they all lie on the border
#        while count_limit < 100:
##           or the code gets stuck because of isolated regions
#            r_adjoin, c_adjoin = adjoin.nonzero()
##           Obtain coordinates of adjoining unwrapped phase pixels
#            if r_adjoin.shape == adjoin_stuck:
#                count_limit += 1
##               Make sure loop doesn't get stuck
#            else:
#                count_limit = 0
#            temp = r_adjoin.shape
#            adjoin_stuck = temp
#            for i in range(1, temp[1] + 1):
#                r_active = r_adjoin[i]
#                c_active = c_adjoin[i]
#                if (r_active <= r_dim - 1 and r_active >= 2 and c_active <= c_dim - 1
#                    and c_active >= 2):
##                   Ignore pixels near the border
##                   First search below for an adjoining unwrapped phase pixel
#                    if (branch_cuts[r_active + 1, c_active] == 0
#                        and unwrapped_binary[r_active + 1, c_active] == 1):
#                        phase_ref = IM_unwrapped[r_active + 1, c_active]
##                       Obtain the reference unwrapped phase
#                        p = scipy.unwrap([phase_ref, IM_phase(r_active, c_active)])
##                       Unwrap the active pixel
#                        IM_unwrapped[r_active, c_active] = p[1]
#                        unwrapped_binary[r_active, c_active] = 1
##                       Mark the pixel as unwrapped
#                        adjoin(r_active, c_active) = 0
##                       Remove it from the list of adjoining pixels
#                        if (r_active - 1 >= 1
#                            and unwrapped_binary[r_active - 1, c_active] == 0
#                            and branch_cuts[r_active - 1, c_active] == 0):
#                            adjoin[r_active - 1, c_active] = 1
#                        if (c_active - 1 >= 1
#                            and unwrapped_binary[r_active, c_active - 1] == 0
#                            and branch_cuts[r_active, c_active - 1] == 0):
#                            adjoin[r_active, c_active - 1] = 1
#                        if (c_active + 1 <= c_dim
#                            and unwrapped_binary[r_active, c_active + 1] == 0
#                            and branch_cuts[r_active, c_active + 1] == 0):
#                            adjoin[r_active, c_active + 1] = 1
##                   Then search above
#                    if (branch_cuts[r_active - 1, c_active] == 0
#                        and unwrapped_binary[r_active-1, c_active] == 1):
#                        phase_ref = IM_unwrapped[r_active - 1, c_active]
##                       Obtain the reference unwrapped phase
#                        p = unwrap([phase_ref, IM_phase[r_active, c_active]])
##                       Unwrap the active pixel
#                        IM_unwrapped[r_active, c_active] = p[1]
#                        unwrapped_binary[r_active,c_active] = 1
##                       Mark the pixel as unwrapped
#                        adjoin[r_active, c_active] = 0
##                       Remove it from the list of adjoining pixels
##                       Update the new adjoining pixels:
#                        if (r_active + 1 <= r_dim
#                            and unwrapped_binary[r_active + 1, c_active] == 0
#                            and branch_cuts[r_active + 1, c_active] == 0):
#                            adjoin(r_active + 1, c_active) = 1
#                        if (c_active - 1 >= 1
#                            and unwrapped_binary[r_active, c_active - 1] == 0
#                            and branch_cuts[r_active, c_active - 1] == 0):
#                            adjoin[r_active, c_active - 1] = 1
#                        if (c_active + 1 <= c_dim 
#                            and unwrapped_binary[r_active, c_active + 1] == 0
#                            and branch_cuts[r_active, c_active + 1] == 0):
#                            adjoin[r_active, c_active + 1] = 1
##                   Then search on the right
#                    if (branch_cuts[r_active, c_active + 1] == 0
#                        and unwrapped_binary[r_active, c_active + 1] == 1):
#                        phase_ref = IM_unwrapped[r_active, c_active + 1]
##                       Obtain the reference unwrapped phase
#                        p = unwrap([phase_ref, IM_phase[r_active, c_active]])
##                       Unwrap the active pixel
#                        IM_unwrapped[r_active, c_active] = p[1]
#                        unwrapped_binary[r_active, c_active] = 1
##                       Mark the pixel as unwrapped
#                        adjoin(r_active, c_active) = 0
##                       Remove it from the list of adjoining pixels
#                        if (r_active + 1 <= r_dim
#                            and unwrapped_binary[r_active + 1, c_active] == 0
#                            and branch_cuts[r_active + 1, c_active] == 0):
#                            adjoin[r_active + 1, c_active] = 1
#                        if (c_active - 1 >= 1
#                            and unwrapped_binary[r_active, c_active - 1] == 0
#                            and branch_cuts[r_active, c_active - 1] == 0):
#                            adjoin[r_active, c_active - 1] = 1
#                        if (r_active - 1 >= 1 
#                            and unwrapped_binary[r_active - 1, c_active] == 0
#                            and branch_cuts[r_active - 1, c_active] == 0):
#                            adjoin[r_active - 1, c_active] = 1
##                   Finally search on the left
#                    if (branch_cuts[r_active, c_active - 1] == 0
#                        and unwrapped_binary[r_active, c_active - 1] == 1):
#                        phase_ref = IM_unwrapped[r_active, c_active - 1]
##                       Obtain the reference unwrapped phase
#                        p = unwrap([phase_ref, IM_phase[r_active, c_active]])
##                       Unwrap the active pixel
#                        IM_unwrapped[r_active, c_active] = p(1)
#                        unwrapped_binary(r_active, c_active) = 1
##                       Mark the pixel as unwrapped
#                        adjoin[r_active, c_active] = 0
##                       Remove it from the list of adjoining pixels
#                        if (r_active + 1 <= r_dim
#                            and unwrapped_binary[r_active + 1, c_active] == 0
#                            and branch_cuts[r_active + 1, c_active] == 0):
#                            adjoin[r_active + 1, c_active] = 1
#                        if (c_active + 1 <= c_dim
#                            and unwrapped_binary[r_active, c_active + 1] == 0
#                            and branch_cuts[r_active, c_active + 1] == 0):
#                            adjoin[r_active, c_active + 1] = 1
#                        if (r_active - 1 >= 1
#                            and unwrapped_binary[r_active - 1, c_active] == 0
#                            and branch_cuts[r_active - 1, c_active] == 0):
#                            adjoin[r_active - 1, c_active] = 1

## Finally, fill in the branch cut pixels that adjoin the unwrapped pixels.
## This can be done because the branch cuts actually lie between the pixels,
## and not on top of them.
#    print('Filling in branch cuts that border on unwrapped pixels ...')
#    adjoin = np.zeros((r_dim, c_dim))
##   Re-load the adjoining pixel matrix with the branch cut values:
#    for i in range(2, r_dim):
#        for j in range(2, c_dim):
#           # Identify which branch cut pixel borders an unwrapped pixel
#           if (branch_cuts[i, j] == 1
#              and (
#                branch_cuts[i + 1,j] == 0
#                or branch_cuts[i - 1, j] == 0
#                or branch_cuts[i, j - 1] == 0
#                or branch_cuts[i, j + 1] == 0)):
#             adjoin[i, j] = 1

#    [r_adjoin, c_adjoin] = find(adjoin)
##   Obtain coordinates of adjoining unwrapped phase pixels
#    temp=size(r_adjoin);
#    for i in range(1, temp[1] + 1):
#        r_active = r_adjoin[i]
#        c_active = c_adjoin[i]
##       First search below for an adjoining unwrapped phase pixel
#        if unwrapped_binary(r_active+1, c_active)==1
#            phase_ref=IM_unwrapped(r_active+1, c_active);                                   %Obtain the reference unwrapped phase
#            p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#            IM_unwrapped(r_active, c_active)=p(2);
#            unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#            adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        end
#        %Then search above
#        if unwrapped_binary(r_active-1, c_active)==1
#            phase_ref=IM_unwrapped(r_active-1, c_active);                                   %Obtain the reference unwrapped phase
#            p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#            IM_unwrapped(r_active, c_active)=p(2);
#            unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#            adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        end
#        %Then search on the right
#        if unwrapped_binary(r_active, c_active+1)==1
#            phase_ref=IM_unwrapped(r_active, c_active+1);                                   %Obtain the reference unwrapped phase
#            p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#            IM_unwrapped(r_active, c_active)=p(2);
#            unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#            adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        end
#        %Finally search on the left
#        if unwrapped_binary(r_active, c_active-1)==1
#            phase_ref=IM_unwrapped(r_active, c_active-1);                                   %Obtain the reference unwrapped phase
#            p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#            IM_unwrapped(r_active, c_active)=p(2);
#            unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#            adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        end
#end
#%figure; imagesc(adjoin), colormap(gray), axis square, axis off, title('Branch cut adjoining pixels');
#%figure; imagesc(IM_unwrapped), colormap(gray), axis square, axis off, title('Peripheral branch cut pixels unwrapped');

#t=toc;
#disp(['Floodfill operation completed in ',int2str(t),' second(s).']);
#disp('Done!');



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% GoldsteinUnwrap2D implements 2D Goldstein branch cut phase unwrapping algorithm.
#%
#% References::
#% 1. R. M. Goldstein, H. A. Zebken, and C. L. Werner, Satellite radar interferometry:
#%    Two-dimensional phase unwrapping, Radio Sci., vol. 23, no. 4, pp. 713720, 1988.
#% 2. D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
#%    Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
#%
#% Inputs: 1. Complex image in .mat double format
#%         2. Binary mask (optional)
#% Outputs: 1. Unwrapped phase image
#%          2. Phase quality map
#%
#% This code can easily be extended for 3D phase unwrapping.
#% Posted by Bruce Spottiswoode on 22 December 2008
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% REPLACE WITH YOUR IMAGES
#load 'IM.mat'                               %Load complex image
#IM_mask=ones(size(IM));                     %Mask (if applicable)
#%%
#IM_mag=abs(IM);                             %Magnitude image
#IM_phase=angle(IM);                         %Phase image

#%%  Set parameters
#max_box_radius=4;                           %Maximum search box radius (pixels)
#threshold_std=5;                            %Number of noise standard deviations used for thresholding the magnitude image

#%% Unwrap
#residue_charge=PhaseResidues(IM_phase, IM_mask);                            %Calculate phase residues
#branch_cuts=BranchCuts(residue_charge, max_box_radius, IM_mask);            %Place branch cuts
#[IM_unwrapped, rowref, colref]=FloodFill(IM_phase, branch_cuts, IM_mask);   %Flood fill phase unwrapping

#%% Display results
#figure; imagesc(residue_charge), colormap(gray), axis square, axis off, title('Phase residues (charged)');
#figure; imagesc(branch_cuts), colormap(gray), axis square, axis off, title('Branch cuts');
#figure; imagesc(immultiply(IM_phase,IM_mask)), colormap(gray), axis square, axis off, title('Wrapped phase');
#tempmin=min(min(IM_unwrapped));          %This bit is just done to create a pleasing display when a mask is used
#temp=(IM_unwrapped==0);
#temp_IM=IM_unwrapped;
#temp_IM(temp)=tempmin;
#figure; imagesc(temp_IM), colormap(gray), axis square, axis off, title('Unwrapped phase');



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% GuidedFloodFill.m unwraps a single 2D image.
#% Input: IM_phase, IM_unwrapped (seed points / pixels already unwrapped),
#% unwrapped_binary the derivative variance, an adjoining matrix and a mask.
#% Created by B.S. Spottiswoode on 11/11/2004
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#function IM_unwrapped=GuidedFloodFill(IM_phase, IM_unwrapped, unwrapped_binary, derivative_variance, adjoin, IM_mask)

#[r_dim, c_dim]=size(IM_phase);
#isolated_adjoining_pixel_flag=1;                            %Only remains set if an isolated adjoining pixel exists

#while sum(sum(adjoin(2:r_dim-1,2:c_dim-1)))~=0              %Loop until there are no more adjoining pixels
#adjoining_derivative_variance=derivative_variance.*adjoin + 100.*~adjoin;      %Derivative variance values of the adjoining pixels (pad the zero adjoining values with 100)
#[r_adjoin, c_adjoin]=find(adjoining_derivative_variance==min(min(adjoining_derivative_variance)));      %Obtain coordinates of the maximum adjoining unwrapped phase pixel
#r_active=r_adjoin(1);
#c_active=c_adjoin(1);
#isolated_adjoining_pixel_flag=1;                                                %This gets cleared as soon as the pixel is unwrapped
#if (r_active==r_dim | r_active==1 | c_active==c_dim | c_active==1);             %Ignore pixels near the border
#   IM_unwrapped(r_active, c_active)=0;
#   adjoin(r_active, c_active)=0;
#else
#    %First search below for an adjoining unwrapped phase pixel
#    if unwrapped_binary(r_active+1, c_active)==1
#        phase_ref=IM_unwrapped(r_active+1, c_active);                                   %Obtain the reference unwrapped phase
#        p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#        IM_unwrapped(r_active, c_active)=p(2);
#        unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#        adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        isolated_adjoining_pixel_flag=0;
#        %Update the new adjoining pixels:
#        if r_active-1>=1 & unwrapped_binary(r_active-1, c_active)==0 & IM_mask(r_active-1, c_active)==1
#            adjoin(r_active-1, c_active)=1;
#        end
#        if c_active-1>=1 & unwrapped_binary(r_active, c_active-1)==0 & IM_mask(r_active, c_active-1)==1
#            adjoin(r_active, c_active-1)=1;
#        end
#        if c_active+1<=c_dim & unwrapped_binary(r_active, c_active+1)==0 & IM_mask(r_active, c_active+1)==1
#            adjoin(r_active, c_active+1)=1;
#        end
#    end
#    %Then search above
#    if unwrapped_binary(r_active-1, c_active)==1
#        phase_ref=IM_unwrapped(r_active-1, c_active);                                       %Obtain the reference unwrapped phase
#        p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#        IM_unwrapped(r_active, c_active)=p(2);
#        unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#        adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        isolated_adjoining_pixel_flag=0;
#        %Update the new adjoining pixels:
#        if r_active+1<=r_dim & unwrapped_binary(r_active+1, c_active)==0 & IM_mask(r_active+1, c_active)==1
#            adjoin(r_active+1, c_active)=1;
#        end
#        if c_active-1>=1 & unwrapped_binary(r_active, c_active-1)==0 & IM_mask(r_active, c_active-1)==1
#            adjoin(r_active, c_active-1)=1;
#        end
#        if c_active+1<=c_dim & unwrapped_binary(r_active, c_active+1)==0 & IM_mask(r_active, c_active+1)==1
#            adjoin(r_active, c_active+1)=1;
#        end
#    end
#    %Then search on the right
#    if unwrapped_binary(r_active, c_active+1)==1
#        phase_ref=IM_unwrapped(r_active, c_active+1);                                       %Obtain the reference unwrapped phase
#        p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#        IM_unwrapped(r_active, c_active)=p(2);
#        unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#        adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        isolated_adjoining_pixel_flag=0;
#        %Update the new adjoining pixels:
#        if r_active+1<=r_dim & unwrapped_binary(r_active+1, c_active)==0 & IM_mask(r_active+1, c_active)==1
#            adjoin(r_active+1, c_active)=1;
#        end
#        if c_active-1>=1 & unwrapped_binary(r_active, c_active-1)==0 & IM_mask(r_active, c_active-1)==1
#            adjoin(r_active, c_active-1)=1;
#        end
#        if r_active-1>=1 & unwrapped_binary(r_active-1, c_active)==0 & IM_mask(r_active-1, c_active)==1
#            adjoin(r_active-1, c_active)=1;
#        end
#    end
#    %Finally search on the left
#    if unwrapped_binary(r_active, c_active-1)==1
#        phase_ref=IM_unwrapped(r_active, c_active-1);                                       %Obtain the reference unwrapped phase
#        p=unwrap([phase_ref IM_phase(r_active, c_active)]);                             %Unwrap the active pixel
#        IM_unwrapped(r_active, c_active)=p(2);
#        unwrapped_binary(r_active, c_active)=1;                                         %Mark the pixel as unwrapped
#        adjoin(r_active, c_active)=0;                                                   %Remove it from the list of adjoining pixels
#        isolated_adjoining_pixel_flag=0;
#        %Update the new adjoining pixels:
#        if r_active+1<=r_dim & unwrapped_binary(r_active+1, c_active)==0 & IM_mask(r_active+1, c_active)==1
#            adjoin(r_active+1, c_active)=1;
#        end
#        if c_active+1<=c_dim & unwrapped_binary(r_active, c_active+1)==0 & IM_mask(r_active, c_active+1)==1
#            adjoin(r_active, c_active+1)=1;
#        end
#        if r_active-1>=1 & unwrapped_binary(r_active-1, c_active)==0 & IM_mask(r_active-1, c_active)==1
#            adjoin(r_active-1, c_active)=1;
#        end
#    end
#    if isolated_adjoining_pixel_flag==1;
#        adjoin(r_active,c_active)=0;                                                    %Remove the current active pixel from the adjoin list
#    end
#end
#end


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% PhaseDerivativeVariance creates a phase quality map by
#% computing the variance of the partial derivatives of the
#% locally unwrapped phase. This is then used
#% to guide the phase unwrapping path. Uses only the 4 nearest
#% neighbours. The user may also input a binary mask.
#% Created by B.S. Spottiswoode on 18/10/2004
#% Last modified on 06/12/2004
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#function derivative_variance=PhaseDerivativeVariance(IM_phase, varargin);

#[r_dim,c_dim]=size(IM_phase);
#if nargin>=2                                    %Has a mask been included? If so crop the image to the mask borders to save computational time
#    IM_mask=varargin{1};
#    [maskrows,maskcols]=find(IM_mask);          %Identify coordinates of the mask
#    minrow=min(maskrows)-1;                     %Identify the limits of the mask
#    maxrow=max(maskrows)+1;
#    mincol=min(maskcols)-1;
#    maxcol=max(maskcols)+1;
#    width=maxcol-mincol;                        %Now ensure that the cropped area is square
#    height=maxrow-minrow;
#    if height>width
#        maxcol=maxcol + floor((height-width)/2) + mod(height-width,2);
#        mincol=mincol - floor((height-width)/2);
#    elseif width>height
#        maxrow=maxrow + floor((width-height)/2) + mod(width-height,2);
#        minrow=minrow - floor((width-height)/2);
#    end
#    if minrow<1 minrow=1; end
#    if maxrow>r_dim maxrow=r_dim; end
#    if mincol<1 mincol=1; end
#    if maxcol>c_dim maxcol=c_dim; end
#    IM_phase=IM_phase(minrow:maxrow, mincol:maxcol);        %Crop the original image to save computation time
#end

#[dimx, dimy]=size(IM_phase);
#dx=zeros(dimx,dimy);
#p=unwrap([IM_phase(:,1) IM_phase(:,2)],[],2);
#dx(:,1)=(p(:,2) - IM_phase(:,1))./2;                    %Take the partial derivative of the unwrapped phase in the x-direction for the first column
#p=unwrap([IM_phase(:,dimy-1) IM_phase(:,dimy)],[],2);
#dx(:,dimy)=(p(:,2) - IM_phase(:,dimy-1))./2;            %Take the partial derivative of the unwrapped phase in the x-direction for the last column
#for i=2:dimy-1
#    p=unwrap([IM_phase(:,i-1) IM_phase(:,i+1)],[],2);
#    dx(:,i)=(p(:,2) - IM_phase(:,i-1))./3;              %Take partial derivative of the unwrapped phase in the x-direction for the remaining columns
#end

#dy=zeros(dimx,dimy);
#q=unwrap([IM_phase(1,:)' IM_phase(2,:)'],[],2);
#dy(1,:)=(q(:,2)' - IM_phase(1,:))./2;                   %Take the partial derivative of the unwrapped phase in the y-direction for the first row
#p=unwrap([IM_phase(dimx-1,:)' IM_phase(dimx,:)'],[],2);
#dy(dimx,:)=(q(:,2)' - IM_phase(dimx-1,:))./2;           %Take the partial derivative of the unwrapped phase in the y-direction for the last row
#for i=2:dimx-1
#    q=unwrap([IM_phase(i-1,:)' IM_phase(i+1,:)'],[],2);
#    dy(i,:)=(q(:,2)' - IM_phase(i-1,:))./3;             %Take the partial derivative of the unwrapped phase in the y-direction for the remaining rows
#end

#dx_centre=dx(2:dimx-1, 2:dimy-1);
#dx_left=dx(2:dimx-1,1:dimy-2);
#dx_right=dx(2:dimx-1,3:dimy);
#dx_above=dx(1:dimx-2,2:dimy-1);
#dx_below=dx(3:dimx,2:dimy-1);
#mean_dx=(dx_centre+dx_left+dx_right+dx_above+dx_below)./5;

#dy_centre=dy(2:dimx-1, 2:dimy-1);
#dy_left=dy(2:dimx-1,1:dimy-2);
#dy_right=dy(2:dimx-1,3:dimy);
#dy_above=dy(1:dimx-2,2:dimy-1);
#dy_below=dy(3:dimx,2:dimy-1);
#mean_dy=(dy_centre+dy_left+dy_right+dy_above+dy_below)./5;

#stdvarx=sqrt( (dx_left - mean_dx).^2 + (dx_right - mean_dx).^2 + ...
#              (dx_above - mean_dx).^2 + (dx_below - mean_dx).^2 + (dx_centre - mean_dx).^2 );
#stdvary=sqrt( (dy_left - mean_dy).^2 + (dy_right - mean_dy).^2 + ...
#              (dy_above - mean_dy).^2 + (dy_below - mean_dy).^2 + (dy_centre - mean_dy).^2 );
#derivative_variance=100*ones(dimx, dimy);                         %Ensure that the border pixels have high derivative variance values
#derivative_variance(2:dimx-1, 2:dimy-1)=stdvarx + stdvary;

#if nargin>=2                                                      %Does the image have to be padded back to the original size?
#    [orig_rows, orig_cols]=size(IM_mask);
#    temp=100*ones(orig_rows, orig_cols);
#    temp(minrow:maxrow, mincol:maxcol)=derivative_variance;       %Pad the remaining pixels with poor phase quality values
#    derivative_variance=temp;
#end

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% PhaseResidues.m calculates the phase residues for a given wrapped phase
#% image. Note that by convention the positions of the phase residues are
#% marked on the top left corner of the 2 by 2 regions.
#%
#%   active---res4---right
#%      |              |
#%     res1           res3
#%      |              |
#%   below---res2---belowright
#% Phase residues with integer multiples of 2*pi are not accounted for, but
#% these rarely occur.
#% Created by B.S. Spottiswoode on 07/10/2004
#% Last modified on 08/10/2004
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#function residue_charge=PhaseResidues(IM_phase, IM_mask);

#[rows, cols]=size(IM_phase);

#%The code below is simply a vectorised representation of operations on 2 by 2
#%blocks in the matrix
#IM_active=IM_phase;
#IM_below=zeros(rows,cols);
#IM_below(1:rows-1,:)=IM_phase(2:rows,:);
#IM_right=zeros(rows,cols);
#IM_right(:,1:cols-1)=IM_phase(:,2:cols);
#IM_belowright=zeros(rows,cols);
#IM_belowright(1:rows-1,1:cols-1)=IM_phase(2:rows,2:cols);

#res1=mod(IM_active - IM_below + pi, 2*pi) - pi;          %Wrap the phase differences as we loop around the 2 by 2 blocks
#res2=mod(IM_below - IM_belowright + pi, 2*pi) - pi;
#res3=mod(IM_belowright - IM_right + pi, 2*pi) - pi;
#res4=mod(IM_right - IM_active + pi, 2*pi) - pi;

#temp_residues=res1+res2+res3+res4;              %Sum the phase differences. Positive residues appear as 2*pi, negative as -2*pi.
#residues=(temp_residues>=6);                    %Assign 1 to positive residue (which should equal 2*pi)
#residues=residues - (temp_residues<=-6);        %Assign -1 to negative residues (which should equal -2*pi)
#residues(:,cols)=0; residues(rows,:)=0;         %Zero pad the border residues
#residues(:,1)=0; residues(1,:)=0;
#residue_charge=residues;

#residue_sum=sum(sum(abs(residues)));




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% QualityGuidedUnwrap2D implements 2D quality guided path following phase
#% unwrapping algorithm.
#%
#% Technique adapted from:
#% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
#% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
#%
#% Inputs: 1. Complex image in .mat double format
#%         2. Binary mask (optional)
#% Outputs: 1. Unwrapped phase image
#%          2. Phase quality map
#%
#% This code can easily be extended for 3D phase unwrapping.
#% Posted by Bruce Spottiswoode on 22 December 2008
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% REPLACE WITH YOUR IMAGES
#load 'IM.mat'                               %Load complex image
#im_mask=ones(size(IM));                     %Mask (if applicable)
#%%

#im_mag=abs(IM);                             %Magnitude image
#im_phase=angle(IM);                         %Phase image
#im_unwrapped=zeros(size(IM));               %Zero starting matrix for unwrapped phase
#adjoin=zeros(size(IM));                     %Zero starting matrix for adjoin matrix
#unwrapped_binary=zeros(size(IM));           %Binary image to mark unwrapped pixels

#%% Calculate phase quality map
#im_phase_quality=PhaseDerivativeVariance(im_phase);

#%% Identify starting seed point on a phase quality map
#minp=im_phase_quality(2:end-1, 2:end-1); minp=min(minp(:));
#maxp=im_phase_quality(2:end-1, 2:end-1); maxp=max(maxp(:));
#figure; imagesc(im_phase_quality,[minp maxp]), colormap(gray), axis square, axis off; title('Phase quality map');
#uiwait(msgbox('Select known true phase reference phase point. Black = high quality phase; white = low quality phase.','Phase reference point','modal'));
#[xpoint,ypoint] = ginput(1);                %Select starting point for the guided floodfill algorithm

#%% Unwrap
#colref=round(xpoint); rowref=round(ypoint);
#im_unwrapped(rowref,colref)=im_phase(rowref,colref);                        %Save the unwrapped values
#unwrapped_binary(rowref,colref,1)=1;
#if im_mask(rowref-1, colref, 1)==1 adjoin(rowref-1, colref, 1)=1; end       %Mark the pixels adjoining the selected point
#if im_mask(rowref+1, colref, 1)==1 adjoin(rowref+1, colref, 1)=1; end
#if im_mask(rowref, colref-1, 1)==1 adjoin(rowref, colref-1, 1)=1; end
#if im_mask(rowref, colref+1, 1)==1 adjoin(rowref, colref+1, 1)=1; end
#im_unwrapped=GuidedFloodFill(im_phase, im_unwrapped, unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

#figure; imagesc(im_mag), colormap(gray), axis square, axis off; title('Magnitude image');
#figure; imagesc(im_phase), colormap(gray), axis square, axis off; title('Wrapped phase');
#figure; imagesc(im_unwrapped), colormap(gray), axis square, axis off; title('Unwrapped phase');
#"""
