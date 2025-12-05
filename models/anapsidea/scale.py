import numpy as np

vs = []
ls = []
with open("curve.obj", 'r') as f:
    lines = f.readlines()

    for line in lines:
        if line.startswith("v "):
            rowvs = line.split(" ")[1:4]
            rowvs = [float(v) for v in rowvs]
            vs.append(rowvs)
        elif line.startswith("l "):
            rowls = line.split(" ")[1:]
            rowls = [int(l) for l in rowls]
            ls.append(rowls)

vs = np.array(vs)
ls = np.array(ls)

vs *= 0.5

with open("curve.obj", 'w') as f:
    for v in vs:
        f.write(f"v {v[0]} {v[1]} {v[2]}\n")
    for l in ls:
        string = "l "
        for i in range(len(l)):
            string += f"{l[i]} "
        string += "\n"
        f.write(string)
