#-*- coding: utf-8 -*-
import os
import math
import random
import numpy as np
import subprocess
import sys
import shutil
from shutil import copy


# read ligand coordinate
def read(path, ligtype):
    file = open(path)
    lines = file.readlines()
    act = [[]]
    n = 0
    for i in range(len(lines)):
        line = lines[i]
        if line[:6] == 'HETATM' and line[17:20] == ligtype and line[77] != 'H':
            act[n].append(float(line[30:38]))
            act[n].append(float(line[38:46]))
            act[n].append(float(line[46:54]))
            n = n+1
            act.append([])
    del act[-1]
    return act

# Remove duplicate residues
def repeat(test, actz):
    y = 0
    for c in range(len(actz)):
        if test == actz[c]:
            y = 1
            break
    if y == 0:
        actz.append(test)

# calculate centerPoint coordinate
def centerPoint(array):
    x, y, z=[], [], []
    k = len(array)
    for l in range(k):
        x.append(array[l][0])
        y.append(array[l][1])
        z.append(array[l][2])
    center = []
    sum1,sum2,sum3=0,0,0
    for i in range(k):
        sum1 = sum1+float(x[i])
        sum2 = sum2+float(y[i])
        sum3 = sum3+float(z[i])
    center.append(sum1/len(x))
    center.append(sum2/len(y))
    center.append(sum3/len(z))
    return center

def energyclash(coordinate, actw):
    atpatom = ['P', 'O', 'O', 'O',
               'P', 'O', 'O', 'O',
               'P', 'O', 'O', 'O',
               'O', 'C', 'C', 'O', 'C', 'O', 'C', 'O', 'C',
               'N', 'C', 'N', 'C', 'C', 'N', 'N', 'C', 'N', 'C']
    dict = {'C':1.70, 'N':1.55, 'O':1.52, 'S':1.80, 'P':1.80} # VDW radii
    dict2 = {'C-C': 1.41, 'N-N': 1.46, 'O-O': 1.21, 'P-P': 2.21,
             'C-N': 1.34, 'C-O': 1.43, 'C-S': 1.81, 'C-P': 1.87,
             'N-O': 1.44, 'N-P': 1.77, 'O-S': 1.51, 'O-P': 1.45} # covalent bond radii
    # ligand and protein
    exclash = 0
    for p in range(31):
        for m in range(len(coordinate)):
            d = 0
            for r in range(3):
                d = d + (float(coordinate[m][r + 2]) - float(actw[p][r])) ** 2
            d = d ** 0.5
            c1 = coordinate[m][1]
            code = c1 + '-' + atpatom[p]
            labelis = code in dict2.keys()
            if labelis == False:
                code = atpatom[p] + '-' + c1
                labelis = code in dict2.keys()
                if labelis == False:
                    r1 = dict[c1]
                    r2 = dict[atpatom[p]]
                    base = 0.5 * (r1 + r2)
                    if d < base:  # clash
                        exclash = exclash + (base - d) / base
                else:
                    distance = dict2[code]
                    if d < distance:  # clash
                        exclash = exclash + (distance - d) / distance
            else:
                distance = dict2[code]
                if d < distance:  # clash
                    exclash = exclash + (distance - d) / distance

    # internal ligand
    inclash = 0
    suan = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11],
            [12, 13, 14, 15, 16, 17, 18, 19, 20],
            [21, 22, 23, 24, 25, 26, 27, 28, 29, 30]]
    for t in range(4):
        for u in range(len(suan[t])):
            linea = atpatom[suan[t][u]]
            count = t + 1
            for h in range(4 - t):
                for j in range(len(suan[count])):
                    lineb = atpatom[suan[count][j]]
                    dis = 0
                    for k in range(3):
                        dis = dis + (actw[suan[t][u]][k] - actw[suan[count][j]][k]) ** 2
                    dis = dis ** 0.5
                    code1 = linea + '-' + lineb
                    labelis1 = code1 in dict2.keys()
                    if labelis1 == False:
                        code1 = lineb + '-' + linea
                        distance1 = dict2[code1]
                        if dis < distance1:
                            inclash = inclash + (distance1 - dis) / distance1
                    else:
                        distance1 = dict2[code1]
                        if dis < distance1:
                            inclash = inclash + (distance1 - dis) / distance1

                count = count + 1
    energy = exclash + inclash
    return energy

# The residue index required for clash calculation
def clashResidue(residueIndex, length):
    clash = []
    for i in range(len(residueIndex)):
        for j in range(1, 16):
            c = int(residueIndex[i])-j
            if c > 0:
                repeat(c, clash)
        repeat(int(residueIndex[i]), clash)
        for j in range(1,16):
            c = int(residueIndex[i])+j
            if c < length:
                repeat(c, clash)
    return clash

def energy(ATPdock, PATH):
    global energyscore
    os.chdir(ATPdock+'/basefile')
    p = subprocess.Popen('./pythonsh prepare_ligand4.py -l '+ PATH + '/final.pdb'+' -o '+ PATH + '/final.pdbqt', shell=True)
    p.wait()
    result = subprocess.Popen('./vina --score_only --receptor '+ PATH + '/pdb.pdbqt --ligand '+ PATH + '/final.pdbqt', shell=True, stdout=subprocess.PIPE)
    lines = result.stdout.readlines()
    for i in range(-10, -4):
        if lines[i][0:8] == b'Affinity':
            energyscore = float(lines[i][9:17])
            break
    os.remove(PATH + '/final.pdbqt')
    os.remove(PATH + '/final.pdb')
    return energyscore

def completion(w, move, acty, acty1):
    actComplete = [[0] * 3 for t in range(31)]
    if w == 0 or w == 1 or w == 2 or w == 4 or w == 6:
        j = 0
        for i in range(31):
            if i > (move[-1] - 1):
                actComplete[i] = acty[i]
            else:
                actComplete[i] = acty1[j]
                j = j + 1

    elif w == 3 or w == 5:
        j = 0
        for i in range(31):
            if i > (move[-1] - 1) or i == move[-4]:
                actComplete[i] = acty[i]
            else:
                actComplete[i] = acty1[j]
                j = j + 1

    elif w == 7:
        j = 0
        for i in range(31):
            if i > (move[-1] - 1) or i == 0:
                actComplete[i] = acty[i]
            else:
                actComplete[i] = acty1[j]
                j = j + 1
    return actComplete


def rmsd(actBase, act):
    distance = 0
    for l in range(31):
        chazhi = 0
        for d in range(3):
            chazhi = chazhi + (actBase[l][d] - act[l][d]) ** 2
        distance = distance + chazhi
    distance = distance / 31
    distance = round(distance ** 0.5, 2)
    return distance

# save ligand file
def generatepdb(PATH, correct, acty3, name):
    act = acty3[:]
    data = open(PATH + '/' + name + '.pdb', 'w+')
    i = 0
    file = open(correct)
    lines = file.readlines()
    for n in range(len(lines)):
        line = lines[n]
        if line[:6] == 'HETATM' and line[77] != 'H':
            print(line[0:30], end='', file=data)
            for j in range(3):
                coordinate = str(act[i][j])
                if j == 0:
                    if coordinate[-2] == '.':
                        coordinate = str(act[i][j]) + '00'
                    if coordinate[-3] == '.':
                        coordinate = str(act[i][j]) + '0'
                    if coordinate[4] == '.':
                        print(coordinate, end='', file=data)
                    if coordinate[3] == '.':
                        print('', coordinate, end='', file=data)
                    if coordinate[2] == '.':
                        print(' ', coordinate, end='', file=data)
                    if coordinate[1] == '.':
                        print('  ', coordinate, end='', file=data)

                if j == 1 or j == 2:
                    if coordinate[-2] == '.':
                        coordinate = str(act[i][j]) + '00'
                    if coordinate[-3] == '.':
                        coordinate = str(act[i][j]) + '0'
                    if coordinate[4] == '.':
                        print(coordinate, end='', file=data)
                    if coordinate[3] == '.':
                        print('', coordinate, end='', file=data)
                    if coordinate[2] == '.':
                        print(' ', coordinate, end='', file=data)
                    if coordinate[1] == '.':
                        print('  ', coordinate, end='', file=data)
            i = i + 1
            print(' ', line[56:], end='', file=data)
    data.close()

# rotation and translation of ATP
def out(acty, unit):
    rotateActy = acty[:]

    for w in range(8):
        if w == 0:
            move = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
            rotate = [21, 20]
        if w == 1:
            move = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
            rotate = [14, 13]
        if w == 2:
            move = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
            rotate = [13, 12]
        if w == 3:
            move = [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12]
            rotate = [12, 8]
        if w == 4:
            move = [1, 2, 3, 4, 5, 6, 7, 8]
            rotate = [8, 11]
        if w == 5:
            move = [1, 2, 3, 4, 6, 7, 8]
            rotate = [11, 4]
        if w == 6:
            move = [1, 2, 3, 4]
            rotate = [4, 7]
        if w == 7:
            move = [2, 3, 4]
            rotate = [7, 0]

        static = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28,
                29, 30, 31]
        num = [[0] * 3 for i in range(len(move))]
        actstatic = rotateActy[:]
        for t in range(len(move)):
            num[t] = rotateActy[move[t] - 1]
        for t in range(len(move)):
            a = move[t] - 1 - t
            actstatic.pop(a)
            static.pop(a)

        head = rotateActy[rotate[0]]
        tail = rotateActy[rotate[1]]
        dis = 0
        for t in range(3):
            dis = dis + (head[t] - tail[t]) ** 2
        length = math.sqrt(dis)
        unitVector = []
        for t in range(3):
            unitVector.append((head[t] - tail[t]) / length)

        x = unitVector[0]
        y = unitVector[1]
        z = unitVector[2]
        c = math.cos(unit[w])
        s = math.sin(unit[w])

        R = [[(x ** 2) * (1 - c) + c, y * x * (1 - c) + z * s, x * z * (1 - c) - y * s],
             [x * y * (1 - c) - z * s, (y ** 2) * (1 - c) + c, y * z * (1 - c) + x * s],
             [x * z * (1 - c) + y * s, y * z * (1 - c) - x * s, (z ** 2) * (1 - c) + c]]
        temMove = [[0] * 3 for t in range(len(move))]
        for i in range(len(num)):
            for j in range(3):
                temMove[i][j] = num[i][j] - head[j]
        rotateMove = [[0] * 3 for t in range(len(move))]
        for i in range(len(temMove)):
            rotateMove[i] = np.dot(temMove[i], R)
            for j in range(3):
                rotateMove[i][j] += head[j]
        rotateActy = completion(w, move, rotateActy, rotateMove)

    center = centerPoint(rotateActy)
    acty1 = [[0] * 3 for i in range(len(rotateActy))]
    for j in range(11, 14):
        unit[j] = unit[j] % (2 * math.pi)
    array_x = np.array(
        [[1, 0, 0], [0, math.cos(unit[11]), -math.sin(unit[11])], [0, math.sin(unit[11]), math.cos(unit[11])]])
    array_y = np.array(
        [[math.cos(unit[12]), 0, math.sin(unit[12])], [0, 1, 0], [-math.sin(unit[12]), 0, math.cos(unit[12])]])
    array_z = np.array(
        [[math.cos(unit[13]), -math.sin(unit[13]), 0], [math.sin(unit[13]), math.cos(unit[13]), 0], [0, 0, 1]])
    u = np.dot(array_x, array_y)
    u = np.dot(u, array_z)
    for b in range(len(rotateActy)):
        for e in range(3):
            acty1[b][e] = round(
                unit[e + 8] + center[e] + u[e][0] * (rotateActy[b][0] - center[0]) + u[e][1] * (rotateActy[b][1] - center[1]) + u[e][
                    2] * (rotateActy[b][2] - center[2]), 3)
    return acty1

def newindi(individual):
    newindividual = individual[:]
    action = random.randint(0, 9)
    if action < 8:
        spin = individual[action] + random.uniform(-0.025 * (math.pi), 0.025 * (math.pi))
        newindividual[action] = spin
    elif action == 8:
        for g in range(3):
            tran = individual[g+action] + random.uniform(-0.1, 0.1)
            newindividual[g+action] = tran
    else:
        for g in range(3):
            spin = individual[11+g] + random.uniform(-0.025 * (math.pi), 0.025 * (math.pi))
            newindividual[11+g] = spin
    return newindividual

# receptor sequence
def restype(path):
    seq = ''
    file = open(path)
    lines = file.readlines()
    for y in range(len(lines) - 1):
        if lines[y][0:4] == 'ATOM' and lines[y + 1][0:4] == 'ATOM':
            if y == 0:
                seq = seq + allresidue[lines[y][17:20]]
            else:
                if lines[y][17:20] != lines[y + 1][17:20]:
                    seq = seq + allresidue[lines[y + 1][17:20]]
    return seq

def searchtm(searchengine):
    print('search template')
    
    ###########################################################
    # Modified by Jun Hu at 20210728 [BELOW]
    ###########################################################
    # Original 
    #os.chdir(ATPdock + '/PPS-search/java/src')
    #p = subprocess.Popen('java submit ' + workpath + ' ' + str(sequence_identity) + ' ' + search_type + ' ' + searchengine, shell=True)
    p = subprocess.Popen('java -jar PPSsearch.jar ' + workpath + ' ' + str(sequence_identity) + ' ' + search_type + ' ' + searchengine, shell=True)
    ###########################################################
    # Modified by Jun Hu at 20210728 [ABOVE]
    ###########################################################
    
    p.wait()
    print('template is done')
        
allresidue = {'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I',
           'PRO': 'P', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W', 'SER': 'S',
           'THR': 'T', 'CYS': 'C', 'MET': 'M', 'ASN': 'N', 'GLN': 'Q',
           'ASP': 'D', 'GLU': 'E', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H'}

#path of ATPdock program
ATPdock = os.path.split(os.path.realpath(__file__))[0]
baseATP = ATPdock + '/basefile/4zibA_ATP.pdb'

#userpath
dock_path = sys.argv[1]

#read sequence identity
#first line is sequence identity
#second line is ligand type
file = open(dock_path + '/tem.txt')
lines = file.readlines()
sequence_identity = float(lines[0].strip())
search_type = lines[1].strip()

receptor = dock_path + '/pdb.pdb'  # input protein
sequence = restype(receptor)
havesite = 0  # site number

if os.path.exists(dock_path + '/pdb.site') == 0:
    havesite = 0
else:
    sitefile = dock_path + '/pdb.site'
    file = open(sitefile)
    lines = file.readlines()
    for i in range(len(lines)):
        line = lines[i]
        if len(line) != 0:
            havesite += 1
            try:
                os.mkdir(dock_path + '/ATPa' + str(havesite))
            except:
                shutil.rmtree(dock_path + '/ATPa' + str(havesite))
                os.mkdir(dock_path + '/ATPa' + str(havesite))
            copy(receptor, dock_path + '/ATPa' + str(havesite) + '/pdb.pdb')
            site = open(dock_path + '/ATPa' + str(havesite) + '/pdb.site', 'w+')
            print(line.strip(), end='', file=site)
            site.close()

if havesite != 0:
    rankscore = []
    for s in range(1, havesite + 1):
        workpath = dock_path + '/ATPa' + str(s)
        searchtm('APoc')        
        file = open(workpath + '/search_list.txt')
        lines = file.readlines()
        line = lines[0].split('~')
        name = line[0]
        psscore = float(line[2])
        if os.path.exists(workpath + '/pocket_' + name + '.ali.txt') == 0:
            f = os.listdir(workpath)
            for h in range(len(f)):
                if f[h][0:3] != 'pdb':
                    os.remove(workpath + '/' + f[h])
            searchtm('PPS-align')
            file = open(workpath + '/search_list.txt')
            lines = file.readlines()
            line = lines[0].split('~')
            name = line[0]
            psscore = float(line[2])
        
        rankscore.append(psscore)

        t = [0, 0, 0]
        u = [[0.0] * 3 for u in range(3)]
        file = open(workpath + '/pocket_' + name + '.ali.txt')
        lines = file.readlines()
        r = 0
        for k in range(len(lines)):
            line = lines[k]
            if (line[0].isdigit() and line[1] == ' ') or (line[0] == ' ' and line[1].isdigit()):
                array = []
                j = len(line)
                for m in range(4, j):
                    if line[m - 1] == ' ' and (line[m] == '-' or line[m].isdigit()):
                        array.append(m)
                t[r] = float(line[array[0]:array[1] - 1])
                u[r][0] = float(line[array[1]:array[2] - 1])
                u[r][1] = float(line[array[2]:array[3] - 1])
                u[r][2] = float(line[array[3]:j])
                r = r + 1

        a = np.array(u)
        litype = name[6:9]
        path1 = workpath + '/' + name + '.pdb'
        act = read(path1, litype)

        solve1 = [[0] * 3 for k in range(len(act))]
        for g in range(len(act)):
            a1 = float(act[g][0]) - t[0]
            b1 = float(act[g][1]) - t[1]
            c1 = float(act[g][2]) - t[2]
            b = np.array([[a1], [b1], [c1]])
            v = np.linalg.solve(a, b)
            solve1[g][0] = round(v[0][0], 3)
            solve1[g][1] = round(v[1][0], 3)
            solve1[g][2] = round(v[2][0], 3)
        actc = solve1

        generatepdb(workpath, path1, actc, 'mu')
        p = subprocess.Popen('obabel -ipdb ' + workpath + '/mu.pdb' + ' -omol2 -O ' + workpath + '/mu.mol2', shell=True)
        p.wait()
        os.chdir(ATPdock+'/basefile')
        align = subprocess.Popen('./LSalign ' + ATPdock + '/basefile/4zibA_ATP.mol2 ' + workpath + '/mu.mol2'+' -rf 1 -o '+ workpath + '/initialATP.pdb -a '+ workpath + '/align.txt -m '+ workpath + '/matrix.txt', shell=True, stdout=subprocess.PIPE)
        align.wait()
        lines = align.stdout.readlines()
        output = lines[3].split()
        pcscore = float(output[2])
        rankscore.append(pcscore)

        p = subprocess.Popen('./pythonsh prepare_receptor4.py -r '+ workpath + '/pdb.pdb'+' -o '+ workpath + '/pdb.pdbqt', shell=True)
        p.wait()

        path = workpath + '/initialATP.pdb'  # initial docking ATP
        acty = read(path, 'QUE')
        initialatp = acty[:]
        generatepdb(workpath, baseATP, acty, 'rough')

        # read binding pocket residues
        file = open(workpath + '/pdb.site')
        lines = file.readlines()
        residueIndex = lines[0].strip('\n')
        residueIndex = residueIndex.split(' ')
        for k in range(len(residueIndex)):
            residueIndex[k] = residueIndex[k][1:]

        # divide receptor to pocket and nopocket
        pocketpdb = open(workpath + '/pocket.pdb', 'w+')
        nopocketpdb = open(workpath + '/nopocket.pdb', 'w+')
        # receptor clash residue coordinate
        file = open(receptor)
        lines = file.readlines()
        clashIndex = clashResidue(residueIndex, len(lines))
        coordinate = [[]]
        n = 0
        for k in range(len(lines)):
            line = lines[k]
            if line[0:4] == 'ATOM' and line[22:26].strip() in residueIndex :
                print(line, end='', file = pocketpdb)
            if line[0:4] == 'ATOM' and line[22:26].strip() not in residueIndex:
                print(line, end='', file = nopocketpdb)
            if line[0:4] == 'ATOM' and line[13] != 'E' and int(line[22:26].strip()) in clashIndex:
                coordinate[n].append(int(line[22:26].strip()))  # residue index
                coordinate[n].append(line[13])  # atom type
                coordinate[n].append(line[30:38])
                coordinate[n].append(line[38:46])
                coordinate[n].append(line[46:54])
                n = n + 1
                coordinate.append([])
        del coordinate[-1]
        pocketpdb.close()
        nopocketpdb.close()

        Tem = 0.3
        Temper = 0.02
        storescore = []
        individual = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        # template score
        generatepdb(workpath, baseATP, initialatp, 'final')
        bestresult = acty[:]
        initialclash = energyclash(coordinate, bestresult)

        if initialclash != 0:
            print('begin docking')
            power = energy(ATPdock, workpath)
            storescore.append(power)
            t = 0
            while True:
                for h in range(10):
                    newindividual = newindi(individual)
                    acty2 = out(acty, newindividual)
                    generatepdb(workpath, baseATP, acty2, 'final')
                    power = energy(ATPdock, workpath)
                    storescore.append(power)
                    receive = 0

                    if storescore[1] < storescore[0]:
                        storescore.pop(0)
                        bestresult = acty2[:]
                        individual = newindividual[:]
                        receive = 1
                        break
                    else:
                        accept = math.exp((storescore[0] - storescore[1]) / Tem)
                        r = random.uniform(0, 1)
                        if r < accept:
                            storescore.pop(0)
                            bestresult = acty2[:]
                            individual = newindividual[:]
                            receive = 1
                            break
                        else:
                            storescore.pop(1)

                if receive == 0:
                    storescore.append(power)
                    rm1 = rmsd(acty, acty2)
                    rm2 = rmsd(acty, bestresult)
                    if rm1 < rm2:
                        storescore.pop(0)
                        bestresult = acty2[:]
                        individual = newindividual[:]
                        receive = 1
                    else:
                        acc = math.exp((rm2 - rm1) / Temper)
                        r = random.uniform(0, 1)
                        if r < acc:
                            storescore.pop(0)
                            bestresult = acty2[:]
                            individual = newindividual[:]
                            receive = 1
                        else:
                            storescore.pop(1)
                if receive == 1:
                    t = t + 1
                if t >= 100:
                    initial1 = energyclash(coordinate, bestresult)
                    if initial1 == 0 or t==30000:
                        break
            generatepdb(workpath, baseATP, bestresult, 'final')
        print('docking is done')
        score = open(workpath + '/score.txt', 'w+')        
        print('PSscore=' + str(psscore), file=score)
        print('PCscore=' + str(pcscore), file=score)
        score.close()

    if len(rankscore) == 2:                
        os.rename(dock_path + '/ATPa1', dock_path + '/ATP1')
    else:
        rank = [[]]
        for i in range(int(len(rankscore) / 2)):
            rank[i].append((rankscore[2*i] + rankscore[2*i+1]) / 2)
            rank[i].append(i+1)
            rank.append([])
        del rank[-1]
        rank.sort()

        for i in range(len(rank)):                        
            os.rename(dock_path + '/ATPa' + str(rank[i][1]), dock_path + '/ATP' + str(len(rank) - i))
    print('ATPdock finished successfully')

else:
    print('There have no pdb.site information, program terminate')