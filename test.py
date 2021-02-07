from multiprocessing import Pool

def square(x):
    for idx in range(100000000):
        var = 0.0
        var += idx**0.00412
    return var

def cube(y):
    for idx in range(100000000):
        var = 0.0
        var += idx**0.01412
    return var

pool = Pool(processes=2)

result_squares = pool.map_async(square, range(1))
result_cubes = pool.map_async(cube, range(1))

print (result_squares.get())
print (result_cubes.get())
