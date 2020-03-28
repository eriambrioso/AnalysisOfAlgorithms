# Recursive implementation of Binary Search

import random
import math
from timeit import default_timer as timer
import matplotlib.pylab as plt


def binarySearch(arr, left, right, x):
    # base case
    if right >= left:
        mid = left + (right - left) // 2
        # element is at the midpoint
        if arr[mid] == x:
            return mid
        # element is in left sub-array
        elif arr[mid] > x:
            return binarySearch(arr, left, mid - 1, x)
        # element is in right sub-array
        else:
            return binarySearch(arr, mid + 1, right, x)
    # element is not found
    return -1


# Driver function
def main():
    sizes = [i*5 for i in range(1, 5)]
    sizes = [i for i in range(10,10000)]
    #print(sizes)
    timing_vals = dict.fromkeys(sizes, 0)
    for size in sizes:
        arr = [i*3 for i in range(size)]
        size = len(arr)
        x = random.randint(0, 3*(size-1))
        start = timer()
        loc = binarySearch(arr, 0, size - 1, x)
        end = timer()
        """
        print("\nResults on an input size of", size, ':')
        if loc == -1:
            print('Element', x, 'not found in list: ')
            print(arr)
        else:
            print('Element', x, 'found at position', loc, 'in list: ')
            print(arr)
        """   
        time_elapsed = end - start

        if time_elapsed > 0.00002:
            time_elapsed =  timing_vals[size-1]
         
        
        timing_vals[size] = time_elapsed

    #print(timing_vals[900], timing_vals[999])
    #x, y = zip(*timing_vals)
    yvals = [6.10365849e-7*math.log2(y) for y in timing_vals.keys()]

    plt.plot(list(timing_vals.keys()), list(timing_vals.values()))
    plt.plot(list(timing_vals.keys()), yvals)
    plt.title("Binary Search time complexity")
    plt.ylabel("time")
    plt.xlabel("n")
    plt.show()


main()
