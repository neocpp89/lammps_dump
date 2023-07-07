#!/usr/bin/env python3
import sys
import math
import pprint

grain_radius = 0.5

class State:
    WAIT_FOR_ITEM = 0
    IN_TIMESTEP = 1
    IN_NUM_ATOMS = 2
    IN_BOX_BOUNDS = 3
    IN_ATOMS = 4

known_items = {
    'ITEM: TIMESTEP': State.IN_TIMESTEP,
    'ITEM: NUMBER OF ATOMS': State.IN_NUM_ATOMS,
    'ITEM: BOX BOUNDS pp pp ss': State.IN_BOX_BOUNDS,
    'ITEM: ATOMS id type x y z vx vy vz': State.IN_ATOMS,
}

frames = []
state = State.WAIT_FOR_ITEM
fr = {}
with open(sys.argv[1]) as f:
    for line in f:
        ls = line.strip()
        if state == State.WAIT_FOR_ITEM:
            if ls in known_items:
                state = known_items[ls]
        elif state == State.IN_TIMESTEP:
            fr['timestep'] = int(ls)
            state = State.WAIT_FOR_ITEM
        elif state == State.IN_NUM_ATOMS:
            fr['num_atoms'] = int(ls)
            state = State.WAIT_FOR_ITEM
        elif state == State.IN_BOX_BOUNDS:
            if 'bounds' not in fr:
                fr['bounds'] = []
            fr['bounds'].append(list(map(float, ls.split())))
            if len(fr['bounds']) == 3:
                state = State.WAIT_FOR_ITEM
        elif state == State.IN_ATOMS:
            toks = ls.split()
            tt = tuple(map(float, toks))
            a = {
                'id': int(tt[0]),
                'type': int(tt[1]),
                'x': tt[2],
                'y': tt[3],
                'z': tt[4],
                'vx': tt[5],
                'vy': tt[6],
                'vz': tt[7],
            }
            if 'atoms' not in fr:
                fr['atoms'] = {}
            fr['atoms'][a['id']] = a
            if len(fr['atoms']) == fr['num_atoms']:
                frames.append(fr)
                fr = {}
                state = State.WAIT_FOR_ITEM

# print('DONE PARSING')

def v_cap(radius, height):
    return math.pi * (height ** 2) * (3 * radius - height) / 3.0

def v_sphere(radius):
    return math.pi * 4.0 * (radius ** 3) / 3.0

def xarea(radius, zsphere, zslice):
    h = max(radius - math.fabs(zsphere - zslice), 0)
    a = math.sqrt(radius ** 2 - (radius - h) ** 2)
    return math.pi * (a ** 2)

def xarea_over_zslice(atoms, z):
    A = 0
    for a in atoms:
        A += xarea(grain_radius, a['z'], z)
    return A

def v_over_zslice(atoms, z):
    vx = 0
    vy = 0
    vz = 0
    for a in atoms:
        A_atom = xarea(grain_radius, a['z'], z)
        vx += (A_atom * a['vx'])
        vy += (A_atom * a['vy'])
        vz += (A_atom * a['vz'])
    A_slice = xarea_over_zslice(atoms, z)
    if A_slice == 0:
        return (0, 0, 0)
    return (vx / A_slice, vy / A_slice, vz / A_slice)

def dvx2_over_zslice(atoms, z):
    dvx2 = 0
    vx, vy, vz = v_over_zslice(atoms, z)
    for a in atoms:
        dvx2 += (xarea(grain_radius, a['z'], z) * ((a['vx']) - vx) ** 2)
    A_slice = xarea_over_zslice(atoms, z)
    if A_slice == 0:
        return 0
    return dvx2 / A_slice

slice_heights = []
slice_dh = 30.0 * grain_radius / 100.0
for i in range(0, 101):
    slice_heights.append(i * slice_dh)

for fr in frames:
    slices = [None] * len(slice_heights)
    for i, s in enumerate(slices):
        slices[i] = dvx2_over_zslice(fr['atoms'].values(), slice_heights[i])
    print(' '.join(f'{s}' for s in slices))
