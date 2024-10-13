HOW TO COMPILE AND RUN THE PROGRAM

1 install gnuplot for generating plots:
    sudo apt-get install gnuplot

2 compile the code:
    gcc galtonbox.c -o galtonbox -lm

3 run the code specifying the number of balls and levels:
    ./galtonbox <number_of_balls> <number_of_levels>

    for example ./galtonbox 50 10

The program generates plots that compare the empirical data with the binomial and normal distributions and a table showing the MSE between empirical and theoretical models for each bin