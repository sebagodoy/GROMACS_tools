#!/usr/bin/env python3
# GROMACS tool
# Repeats cell in gro format

print(" # GROMACS tool: multiply .gro to crystal configuration ")

## Open - parse
InFileName = ''
while InFileName=='':
    InFileName = input('  > Input .gro file: ')
with open(InFileName, 'r') as f:
    InFile = f.readlines()
In_title = InFile[0]
In_NAtoms = int(InFile[1])
In_Cell = [float(i) for i in InFile[In_NAtoms+2].split()]

In_Atoms = []; In_Count_NRes = 0; In_Count_NAtoms = 0; In_Count_AtomsTypes = []
for k in InFile[2:-1]:
    _ = [int(k[:5]),               #ResNumber, int, 5char
         k[5:10],             #ResName, str, 5char
         k[10:15],            #At Type Name, str, 5char
         int(k[15:20]),            #At number total, int, 5char
         [float(k[20:28]),
          float(k[28:36]),
          float(k[36:44])] #At pos, int, 8/8/8
         ]
    # Vels?
    Bool_Vel = False
    if len(k)>44: # velocities are present
        _.append([float(k[44:52]),
                  float(k[52:60]),
                  float(k[60:68])])
        Bool_Vel = True
    # Counters
    if _[0] > In_Count_NRes: In_Count_NRes += 1
    if _[2] not in In_Count_AtomsTypes: In_Count_AtomsTypes.append(_[2])
    In_Count_NAtoms += 1

    In_Atoms.append(_)

# Report
In_NRes = In_Count_NRes
In_AtomTypes = In_Count_AtomsTypes
if not In_NAtoms == In_Count_NAtoms:
    raise IOError(" Number of atoms reported and read does not coincide! Bye.")
print(f'  > Parsed {InFileName} file, got {In_NAtoms} atoms ({len(In_AtomTypes)} types) in {In_NRes} residual groups.')

# Parse base cell (the .gro definition is weird)
Cell = [[In_Cell[0],    .0,         .0],
        [.0,            In_Cell[1], .0],
        [.0,            .0,         In_Cell[2]]]
if len(In_Cell) > 3:
    # values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y)
    # that is:
    # v1(x)              v2(y)           v3(z)
    #  v1(y)=0  v1(z)=0
    #                     v2(x) v2(z)=0
    #                                       v3(x) v3(y)
    # GROMACS only supports boxes with v1(y)=v1(z)=v2(z)=0.
    Cell[1][0] = In_Cell[5] # v2(x)
    Cell[2][0] = In_Cell[7] # v3(x)
    Cell[2][1] = In_Cell[8] # v3(y)
print(' > Base cell is:')
[print(i) for i in Cell]

#### Multiply Cell
CloneAxis = input(' > Clone axis x, y, z as (def: 1 1 1) : ') or '1 1 1'
Xclone, Yclone, Zclone = [int(i) for i in CloneAxis.split()]

print(f' Clonning {Xclone} x {Yclone} x {Zclone} will give '
      f'{Xclone*Yclone*Zclone*In_NAtoms} atoms in '
      f'{Xclone*Yclone*Zclone*In_NRes} residual groups.')

Out_Atoms = []
OutCoun_Atoms = 0; OutCount_Res = 0
for iX in range(Xclone):
    for iY in range(Yclone):
        for iZ in range(Zclone):
            _CountRes = 0 #residuals read reset in each cloned cell
            for iAt in In_Atoms:
                # displace
                xDisp = iX*Cell[0][0] + iY*Cell[1][0] + iZ*Cell[2][0]
                yDisp = iX*Cell[0][1] + iY*Cell[1][1] + iZ*Cell[2][1]
                zDisp = iX*Cell[0][2] + iY*Cell[1][2] + iZ*Cell[2][2]
                # counters
                OutCoun_Atoms +=1
                if not iAt[0] == _CountRes:
                    _CountRes = iAt[0]
                    OutCount_Res +=1
                # build new
                _ = [OutCount_Res,  #res number
                     iAt[1],        # res name
                     iAt[2],        # atom type
                     OutCoun_Atoms, # atom number
                     [m + n for m,n in zip(iAt[4],[xDisp, yDisp, zDisp])]
                     ]
                if len(iAt) == 6: # velocities?
                    _.append(iAt[5])
                Out_Atoms.append(_)

#### Writting out
OutFile = input(' > Out file name (def:(...)_XxYxZ.gro) : ') or ''.join([InFileName.split('.gro')[0],
                                                                         f'_{Xclone}x{Yclone}x{Zclone}.gro'])
with open(OutFile,'w') as f:
    f.write(In_title)
    f.write(str(OutCoun_Atoms).rjust(5)+'\n')
    for iAt in Out_Atoms:
        f.write(str(iAt[0]).rjust(5)) # res number
        f.write(iAt[1]) # res name (5char already)
        f.write(iAt[2]) # At type (5char already)
        f.write(str(iAt[3]).rjust(5)) # At number
        for iCoord in iAt[4]:
            f.write('{:.3f}'.format(iCoord).rjust(8))
        if len(iAt) == 6:
            for iVel in iAt[5]:
                f.write('{:.3f}'.format(iVel).rjust(8))
        f.write('\n')
    # box in nm
    OutBox = []
    for iVect,Mult in zip(Cell, [Xclone, Yclone, Zclone]):
        [OutBox.append('{:.5f}'.format(kUVect*Mult).rjust(10)) for kUVect in iVect]
    f.write(OutBox[0]+OutBox[4]+OutBox[8])
    f.write(OutBox[1]+OutBox[2])
    f.write(OutBox[3]+OutBox[5])
    f.write(OutBox[6]+OutBox[7])
    f.write('\n')

