#from pypower import case9, ppoption, runpf, printpf, runopf, rundcopf
from pypower.api import case9, ppoption, runpf, printpf, runopf, rundcopf
# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

from numpy import \
    ones, zeros, r_, sort, exp, pi, diff, arange, min, \
    argmin, argmax, logical_or, real, imag, any


from sys import stdout, stderr

from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF
from pypower.idx_gen import GEN_BUS, PG, QG, QMAX, QMIN, GEN_STATUS, \
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN
from pypower.idx_brch import F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, \
TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST

from numpy import array
from numpy import flatnonzero as find

from pypower.isload import isload
from pypower.run_userfcn import run_userfcn
from pypower.ppoption import ppoption

def mycase():

    ppc = {"version": '2'}

    ##-----  Power Flow Data  -----##
    ## system MVA base
    ppc["baseMVA"] = 100.0

    ## bus data
    # bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin
    ppc["bus"] = array([
        [1,  3,  0,    0,   0, 0,  1, 1.06,    0,    0, 1, 1.06, 0.94],
        [2,  2, 0, 0, 0, 0,  1, 1.045,  -4.98, 0, 1, 1.06, 0.94],
        [3,  2, 0, 0,   0, 0,  1, 1.01,  -12.72, 0, 1, 1.06, 0.94],
        [4,  1, 0, 0, 0, 0,  1, 1.019, -10.33, 0, 1, 1.06, 0.94],
        [5,  1,  0,  0, 0, 0,  1, 1.02,   -8.78, 0, 1, 1.06, 0.94],
        [6,  2, 0,  0, 0, 0,  1, 1.07,  -14.22, 0, 1, 1.06, 0.94],
        [7,  1,  0,    0,   0, 0,  1, 1.062, -13.37, 0, 1, 1.06, 0.94],
        [8,  2,  0,    0,   0, 0,  1, 1.09,  -13.36, 0, 1, 1.06, 0.94],
        [9,  1, 0, 0, 0, 19, 1, 1.056, -14.94, 0, 1, 1.06, 0.94],
        [10, 1,  0,    0, 0, 0,  1, 1.051, -15.1,  0, 1, 1.06, 0.94],
        [11, 1,  0,  0, 0, 0,  1, 1.057, -14.79, 0, 1, 1.06, 0.94],
        [12, 1,  0,  0, 0, 0,  1, 1.055, -15.07, 0, 1, 1.06, 0.94],
        [13, 1, 0,  0, 0, 0,  1, 1.05,  -15.16, 0, 1, 1.06, 0.94],
        [14, 1, 0,  0,   0, 0,  1, 1.036, -16.04, 0, 1, 1.06, 0.94]
    ])

    ## generator data
    # bus, Pg, Qg, Qmax, Qmin, Vg, mBase, status, Pmax, Pmin, Pc1, Pc2,
    # Qc1min, Qc1max, Qc2min, Qc2max, ramp_agc, ramp_10, ramp_30, ramp_q, apf
    ppc["gen"] = array([
        [1, 232.4, 0, 100,   -100, 1.06,  100, 1, 332.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [2,  40,    0, 50, -100, 1.045, 100, 1, 140,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [3,   0,    0, 40,   -100, 1.01,  100, 1, 100,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [3,   0,    0, 100,   -100, 1,  100, 1, 0,   -100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [4,   0,   0, 100,   -100, 1,  100, 1, 0,   -20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [5,   0,    0, 100,   -100, 1,  100, 1, 0,   -10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [6,   0,    0, 100,  -100, 1.07,  100, 1, 100,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [7,   0,    0, 100,   -100, 1.01,  100, 1, 200,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [8,   0,    0,100,  -100, 1.09,  100, 1, 100,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [9,   0,    0, 100,   -100, 1.01,  100, 1, 0,   -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [10,   0,    0, 100,   -100, 1.01,  100, 1, 0,   -50, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [11,   0,    0, 100,   -100, 1.01,  100, 1, 150,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [12,   0,    0, 100,   -100, 1.01,  100, 1, 0,   -150, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [13,   0,    0, 100,   -100, 1,  100, 1, 0,   -30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [14,   0,    0, 100,   -100, 1.01,  100, 1, 50,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        
    ])

    ## branch data
    # fbus, tbus, r, x, b, rateA, rateB, rateC, ratio, angle, status, angmin, angmax
    ppc["branch"] = array([
        [1,   2, 0.01938, 0.05917, 0.0528, 9900, 0, 0, 0,     0, 1, -360, 360],
        [1,   5, 0.05403, 0.22304, 0.0492, 9900, 0, 0, 0,     0, 1, -360, 360],
        [2,   3, 0.04699, 0.19797, 0.0438, 9900, 0, 0, 0,     0, 1, -360, 360],
        [2,   4, 0.05811, 0.17632, 0.034,  9900, 0, 0, 0,     0, 1, -360, 360],
        [2,   5, 0.05695, 0.17388, 0.0346, 9900, 0, 0, 0,     0, 1, -360, 360],
        [3,   4, 0.06701, 0.17103, 0.0128, 9900, 0, 0, 0,     0, 1, -360, 360],
        [4,   5, 0.01335, 0.04211, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [4,   7, 0,       0.20912, 0,      9900, 0, 0, 0.978, 0, 1, -360, 360],
        [4,   9, 0,       0.55618, 0,      9900, 0, 0, 0.969, 0, 1, -360, 360],
        [5,   6, 0,       0.25202, 0,      9900, 0, 0, 0.932, 0, 1, -360, 360],
        [6,  11, 0.09498, 0.1989,  0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [6,  12, 0.12291, 0.25581, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [6,  13, 0.06615, 0.13027, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [7,   8, 0,       0.17615, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [7,   9, 0,       0.11001, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [9,  10, 0.03181, 0.0845,  0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [9,  14, 0.12711, 0.27038, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [10, 11, 0.08205, 0.19207, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [12, 13, 0.22092, 0.19988, 0,      9900, 0, 0, 0,     0, 1, -360, 360],
        [13, 14, 0.17093, 0.34802, 0,      9900, 0, 0, 0,     0, 1, -360, 360]
    ])

    ##-----  OPF Data  -----##
    ## generator cost data
    # 1 startup shutdown n x1 y1 ... xn yn
    # 2 startup shutdown n c(n-1) ... c0
    ppc["gencost"] = array([
        [2, 0, 0, 3, 0, 20.2, 0],
        [2, 0, 0, 3, 0,      20.1, 0],
        [2, 0, 0, 3, 0,      30.2, 0],
        [2, 0, 0, 3, 0,      40.5, 0],
        [2, 0, 0, 3, 0,      40.6, 0],
        [2, 0, 0, 3, 0, 45.6, 0],
        [2, 0, 0, 3, 0,      21.7, 0],
        [2, 0, 0, 3, 0,      22.8, 0],
        [2, 0, 0, 3, 0,      35.9, 0],
        [2, 0, 0, 3, 0,      50.0, 0],
        [2, 0, 0, 3, 0, 40.2, 0],
        [2, 0, 0, 3, 0,      27.1, 0],
        [2, 0, 0, 3, 0,      43.3, 0],
        [2, 0, 0, 3, 0,      44.4, 0],
        [2, 0, 0, 3, 0,      38.5, 0]
    ])

    return ppc

# Copyright (c) 1996-2015 PSERC. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

"""Prints power flow results.
"""

from sys import stdout

from numpy import \
    ones, zeros, r_, sort, exp, pi, diff, arange, min, \
    argmin, argmax, logical_or, real, imag, any

from numpy import flatnonzero as find

from pypower.idx_bus import BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, \
    VM, VA, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, REF
from pypower.idx_gen import GEN_BUS, PG, QG, QMAX, QMIN, GEN_STATUS, \
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN
from pypower.idx_brch import F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, \
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST

from pypower.isload import isload
from pypower.run_userfcn import run_userfcn
from pypower.ppoption import ppoption


def myprintpf(baseMVA, bus=None, gen=None, branch=None, f=None, success=None,
            et=None, fd=None, ppopt=None):
    """Prints power flow results.

    Prints power flow and optimal power flow results to C{fd} (a file
    descriptor which defaults to C{stdout}), with the details of what
    gets printed controlled by the optional C{ppopt} argument, which is a
    PYPOWER options vector (see L{ppoption} for details).

    The data can either be supplied in a single C{results} dict, or
    in the individual arguments: C{baseMVA}, C{bus}, C{gen}, C{branch}, C{f},
    C{success} and C{et}, where C{f} is the OPF objective function value,
    C{success} is C{True} if the solution converged and C{False} otherwise,
    and C{et} is the elapsed time for the computation in seconds. If C{f} is
    given, it is assumed that the output is from an OPF run, otherwise it is
    assumed to be a simple power flow run.

    Examples::
        ppopt = ppoptions(OUT_GEN=1, OUT_BUS=0, OUT_BRANCH=0)
        fd = open(fname, 'w+b')
        results = runopf(ppc)
        printpf(results)
        printpf(results, fd)
        printpf(results, fd, ppopt)
        printpf(baseMVA, bus, gen, branch, f, success, et)
        printpf(baseMVA, bus, gen, branch, f, success, et, fd)
        printpf(baseMVA, bus, gen, branch, f, success, et, fd, ppopt)
        fd.close()

    @author: Ray Zimmerman (PSERC Cornell)
    """
    ##----- initialization -----
    ## default arguments
    objValue=0
    
    if isinstance(baseMVA, dict):
        have_results_struct = 1
        results = baseMVA
        if gen is None:
            ppopt = ppoption()   ## use default options
        else:
            ppopt = gen
        if (ppopt['OUT_ALL'] == 0):
            return     ## nothin' to see here, bail out now
        if bus is None:
            fd = stdout         ## print to stdout by default
        else:
            fd = bus
        baseMVA, bus, gen, branch, success, et = \
            results["baseMVA"], results["bus"], results["gen"], \
            results["branch"], results["success"], results["et"]
        if 'f' in results:
            f = results["f"]
        else:
            f = None
    else:
        have_results_struct = 0
        if ppopt is None:
            ppopt = ppoption()   ## use default options
            if fd is None:
                fd = stdout         ## print to stdout by default
        if ppopt['OUT_ALL'] == 0:
            return     ## nothin' to see here, bail out now

    isOPF = f is not None    ## FALSE -> only simple PF data, TRUE -> OPF data

    ## options
    isDC            = ppopt['PF_DC']        ## use DC formulation?
    OUT_ALL         = ppopt['OUT_ALL']
    OUT_ANY         = OUT_ALL == 1     ## set to true if any pretty output is to be generated
    OUT_SYS_SUM     = (OUT_ALL == 1) or ((OUT_ALL == -1) and ppopt['OUT_SYS_SUM'])
    OUT_AREA_SUM    = (OUT_ALL == 1) or ((OUT_ALL == -1) and ppopt['OUT_AREA_SUM'])
    OUT_BUS         = (OUT_ALL == 1) or ((OUT_ALL == -1) and ppopt['OUT_BUS'])
    OUT_BRANCH      = (OUT_ALL == 1) or ((OUT_ALL == -1) and ppopt['OUT_BRANCH'])
    OUT_GEN         = (OUT_ALL == 1) or ((OUT_ALL == -1) and ppopt['OUT_GEN'])
    OUT_ANY         = OUT_ANY | ((OUT_ALL == -1) and
                        (OUT_SYS_SUM or OUT_AREA_SUM or OUT_BUS or
                         OUT_BRANCH or OUT_GEN))

    if OUT_ALL == -1:
        OUT_ALL_LIM = ppopt['OUT_ALL_LIM']
    elif OUT_ALL == 1:
        OUT_ALL_LIM = 2
    else:
        OUT_ALL_LIM = 0

    OUT_ANY         = OUT_ANY or (OUT_ALL_LIM >= 1)
    if OUT_ALL_LIM == -1:
        OUT_V_LIM       = ppopt['OUT_V_LIM']
        OUT_LINE_LIM    = ppopt['OUT_LINE_LIM']
        OUT_PG_LIM      = ppopt['OUT_PG_LIM']
        OUT_QG_LIM      = ppopt['OUT_QG_LIM']
    else:
        OUT_V_LIM       = OUT_ALL_LIM
        OUT_LINE_LIM    = OUT_ALL_LIM
        OUT_PG_LIM      = OUT_ALL_LIM
        OUT_QG_LIM      = OUT_ALL_LIM

    OUT_ANY         = OUT_ANY or ((OUT_ALL_LIM == -1) and (OUT_V_LIM or OUT_LINE_LIM or OUT_PG_LIM or OUT_QG_LIM))
    ptol = 1e-4        ## tolerance for displaying shadow prices

    ## create map of external bus numbers to bus indices
    i2e = bus[:, BUS_I].astype(int)
    e2i = zeros(max(i2e) + 1, int)
    e2i[i2e] = arange(bus.shape[0])

    ## sizes of things
    nb = bus.shape[0]      ## number of buses
    nl = branch.shape[0]   ## number of branches
    ng = gen.shape[0]      ## number of generators

    ## zero out some data to make printout consistent for DC case
    if isDC:
        bus[:, r_[QD, BS]]          = zeros((nb, 2))
        gen[:, r_[QG, QMAX, QMIN]]  = zeros((ng, 3))
        branch[:, r_[BR_R, BR_B]]   = zeros((nl, 2))

    ## parameters
    ties = find(bus[e2i[branch[:, F_BUS].astype(int)], BUS_AREA] !=
                   bus[e2i[branch[:, T_BUS].astype(int)], BUS_AREA])
                            ## area inter-ties
    tap = ones(nl)                           ## default tap ratio = 1 for lines
    xfmr = find(branch[:, TAP])           ## indices of transformers
    tap[xfmr] = branch[xfmr, TAP]            ## include transformer tap ratios
    tap = tap * exp(1j * pi / 180 * branch[:, SHIFT]) ## add phase shifters
    nzld = find((bus[:, PD] != 0.0) | (bus[:, QD] != 0.0))
    sorted_areas = sort(bus[:, BUS_AREA])
    ## area numbers
    s_areas = sorted_areas[r_[1, find(diff(sorted_areas)) + 1]]
    nzsh = find((bus[:, GS] != 0.0) | (bus[:, BS] != 0.0))
    allg = find( ~isload(gen) )
    ong  = find( (gen[:, GEN_STATUS] > 0) & ~isload(gen) )
    onld = find( (gen[:, GEN_STATUS] > 0) &  isload(gen) )
    V = bus[:, VM] * exp(-1j * pi / 180 * bus[:, VA])
    out = find(branch[:, BR_STATUS] == 0)        ## out-of-service branches
    nout = len(out)
    if isDC:
        loss = zeros(nl)
    else:
        loss = baseMVA * abs(V[e2i[ branch[:, F_BUS].astype(int) ]] / tap -
                             V[e2i[ branch[:, T_BUS].astype(int) ]])**2 / \
                    (branch[:, BR_R] - 1j * branch[:, BR_X])

    fchg = abs(V[e2i[ branch[:, F_BUS].astype(int) ]] / tap)**2 * branch[:, BR_B] * baseMVA / 2
    tchg = abs(V[e2i[ branch[:, T_BUS].astype(int) ]]      )**2 * branch[:, BR_B] * baseMVA / 2
    loss[out] = zeros(nout)
    fchg[out] = zeros(nout)
    tchg[out] = zeros(nout)

    ##----- print the stuff -----
    if OUT_ANY:
        ## convergence & elapsed time
        if success:
            fd.write('\nConverged in %.2f seconds' % et)
        else:
            fd.write('\nDid not converge (%.2f seconds)\n' % et)

        ## objective function value
        if isOPF:
            fd.write('\nObjective Function Value = %.2f $/hr' % f)

    if OUT_SYS_SUM:
        fd.write('\n================================================================================')
        fd.write('\n|     System Summary                                                           |')
        fd.write('\n================================================================================')
        fd.write('\n\nHow many?                How much?              P (MW)            Q (MVAr)')
        fd.write('\n---------------------    -------------------  -------------  -----------------')
        fd.write('\nBuses         %6d     Total Gen Capacity   %7.1f       %7.1f to %.1f' % (nb, sum(gen[allg, PMAX]), sum(gen[allg, QMIN]), sum(gen[allg, QMAX])))
        fd.write('\nGenerators     %5d     On-line Capacity     %7.1f       %7.1f to %.1f' % (len(allg), sum(gen[ong, PMAX]), sum(gen[ong, QMIN]), sum(gen[ong, QMAX])))
        fd.write('\nCommitted Gens %5d     Generation (actual)  %7.1f           %7.1f' % (len(ong), sum(gen[ong, PG]), sum(gen[ong, QG])))
        fd.write('\nLoads          %5d     Load                 %7.1f           %7.1f' % (len(nzld)+len(onld), sum(bus[nzld, PD])-sum(gen[onld, PG]), sum(bus[nzld, QD])-sum(gen[onld, QG])))
        fd.write('\n  Fixed        %5d       Fixed              %7.1f           %7.1f' % (len(nzld), sum(bus[nzld, PD]), sum(bus[nzld, QD])))
        fd.write('\n  Dispatchable %5d       Dispatchable       %7.1f of %-7.1f%7.1f' % (len(onld), -sum(gen[onld, PG]), -sum(gen[onld, PMIN]), -sum(gen[onld, QG])))
        fd.write('\nShunts         %5d     Shunt (inj)          %7.1f           %7.1f' % (len(nzsh),
            -sum(bus[nzsh, VM]**2 * bus[nzsh, GS]), sum(bus[nzsh, VM]**2 * bus[nzsh, BS]) ))
        fd.write('\nBranches       %5d     Losses (I^2 * Z)     %8.2f          %8.2f' % (nl, sum(loss.real), sum(loss.imag) ))
        fd.write('\nTransformers   %5d     Branch Charging (inj)     -            %7.1f' % (len(xfmr), sum(fchg) + sum(tchg) ))
        fd.write('\nInter-ties     %5d     Total Inter-tie Flow %7.1f           %7.1f' % (len(ties), sum(abs(branch[ties, PF]-branch[ties, PT])) / 2, sum(abs(branch[ties, QF]-branch[ties, QT])) / 2))
        fd.write('\nAreas          %5d' % len(s_areas))
        fd.write('\n')
        fd.write('\n                          Minimum                      Maximum')
        fd.write('\n                 -------------------------  --------------------------------')
        minv = min(bus[:, VM])
        mini = argmin(bus[:, VM])
        maxv = max(bus[:, VM])
        maxi = argmax(bus[:, VM])
        fd.write('\nVoltage Magnitude %7.3f p.u. @ bus %-4d     %7.3f p.u. @ bus %-4d' % (minv, bus[mini, BUS_I], maxv, bus[maxi, BUS_I]))
        minv = min(bus[:, VA])
        mini = argmin(bus[:, VA])
        maxv = max(bus[:, VA])
        maxi = argmax(bus[:, VA])
        fd.write('\nVoltage Angle   %8.2f deg   @ bus %-4d   %8.2f deg   @ bus %-4d' % (minv, bus[mini, BUS_I], maxv, bus[maxi, BUS_I]))
        if not isDC:
            maxv = max(loss.real)
            maxi = argmax(loss.real)
            fd.write('\nP Losses (I^2*R)             -              %8.2f MW    @ line %d-%d' % (maxv, branch[maxi, F_BUS], branch[maxi, T_BUS]))
            maxv = max(loss.imag)
            maxi = argmax(loss.imag)
            fd.write('\nQ Losses (I^2*X)             -              %8.2f MVAr  @ line %d-%d' % (maxv, branch[maxi, F_BUS], branch[maxi, T_BUS]))
        if isOPF:
            minv = min(bus[:, LAM_P])
            mini = argmin(bus[:, LAM_P])
            maxv = max(bus[:, LAM_P])
            maxi = argmax(bus[:, LAM_P])
            fd.write('\nLambda P        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d' % (minv, bus[mini, BUS_I], maxv, bus[maxi, BUS_I]))
            minv = min(bus[:, LAM_Q])
            mini = argmin(bus[:, LAM_Q])
            maxv = max(bus[:, LAM_Q])
            maxi = argmax(bus[:, LAM_Q])
            fd.write('\nLambda Q        %8.2f $/MWh @ bus %-4d   %8.2f $/MWh @ bus %-4d' % (minv, bus[mini, BUS_I], maxv, bus[maxi, BUS_I]))
        fd.write('\n')



    ## bus data
    if OUT_BUS:
        fd.write('\n================================================================================')
        fd.write('\n|     Bus Data                                                                 |')
        fd.write('\n================================================================================')
        fd.write('\n Bus      Voltage          Generation             Load        ')
        if isOPF: fd.write('  Lambda($/MVA-hr)')
        fd.write('\n  #   Mag(pu) Ang(deg)   P (MW)   Q (MVAr)   P (MW)   Q (MVAr)')
        if isOPF: fd.write('     P        Q   ')
        fd.write('\n----- ------- --------  --------  --------  --------  --------')
        if isOPF: fd.write('  -------  -------')
        for i in range(nb):
            fd.write('\n%5d%7.3f%9.3f' % tuple(bus[i, [BUS_I, VM, VA]]))
            if bus[i, BUS_TYPE] == REF:
                fd.write('*')
            else:
                fd.write(' ')
            g  = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) &
                        ~isload(gen))
            ld = find((gen[:, GEN_STATUS] > 0) & (gen[:, GEN_BUS] == bus[i, BUS_I]) &
                        isload(gen))
            if any(g + 1):
                fd.write('%9.2f%10.2f' % (sum(gen[g, PG]), sum(gen[g, QG])))
                objValue = objValue + bus[i, LAM_P] * (sum(gen[g, PG]))
            else:
                fd.write('      -         -  ')

            if logical_or(bus[i, PD], bus[i, QD]) | any(ld + 1):
                if any(ld + 1):
                    fd.write('%10.2f*%9.2f*' % (bus[i, PD] - sum(gen[ld, PG]),
                                                bus[i, QD] - sum(gen[ld, QG])))
                    objValue = objValue + bus[i, LAM_P] * (bus[i, PD] - sum(gen[ld, PG]))
                else:
                    fd.write('%10.2f%10.2f ' % tuple(bus[i, [PD, QD]]))
            else:
                fd.write('       -         -   ')
            if isOPF:
                fd.write('%9.3f' % bus[i, LAM_P])
                if abs(bus[i, LAM_Q]) > ptol:
                    fd.write('%8.3f' % bus[i, LAM_Q])
                else:
                    fd.write('     -')
        fd.write('\n                        --------  --------  --------  --------')
        fd.write('\n               Total: %9.2f %9.2f %9.2f %9.2f' %
            (sum(gen[ong, PG]), sum(gen[ong, QG]),
             sum(bus[nzld, PD]) - sum(gen[onld, PG]),
             sum(bus[nzld, QD]) - sum(gen[onld, QG])))
        fd.write('\n')
        fd.write('===========ObjValue = %f==============\n' %objValue)





def case_iteration():
    ppc = mycase()
    ppopt = ppoption(PF_ALG=2)

    r = runopf(ppc, ppopt)
    
    myprintpf(r)
    
    
    


if __name__ == '__main__':
    case_iteration()
