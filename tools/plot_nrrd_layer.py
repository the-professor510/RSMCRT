import matplotlib.pyplot as plt
import read_nrrd_class

file = "RSMCRT/data/absorb/absorb.nrrd"
read_nrrd = read_nrrd_class.read_nrrd_class()
grid, hdr = read_nrrd.read_nrrd(file)

plt.imshow(grid[100, :, :])
plt.show()

print(sum(sum(grid[100, :, :]))/1000000)