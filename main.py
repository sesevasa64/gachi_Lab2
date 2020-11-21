from figure import *
from mat4d import *


f = Figure('output.txt')
# f += f[1] * Mat4d.parallel(Vec([0, 17, 0]))
# f += f[0] * Mat4d.parallel(Vec([0, 0, 10]))
# f[1].delete(1)
f.plot()
f.save('output.txt')