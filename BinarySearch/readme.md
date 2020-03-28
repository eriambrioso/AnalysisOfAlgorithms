# Binary Search in Python

This is a recursive implementation of binary search in python with some profiling data taken. 

Use the command "python binarySeach.py" to run the program

Below is sample output running on sizes 5, 10, 15, and 20

![Sample output](binsearch_test.PNG)

To collect profiling data I used the timeit function from the timer library and plt() from matplotlib to plot my timing data. The expected time complexity is O(log n), this is graphed in orange. Given the nature of binary search, some iterations took an very long time and caused "super peaks" in the graph. To avoid these peaks skewing the graph, I chose a cap value that helps avoid outliers.

![Timing Data](Figure_1.png)