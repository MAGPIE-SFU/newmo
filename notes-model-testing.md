# Notes

Here are some notes (from CC) about simulating and testing the NFDS model. 

## 2026-09-29

I decided today that the ode solver should give us not the values of the group frequencies
(or counts, if we are using a population of $N$) but instead, the average for the given
year. Since our synthetic carriage population $x(t)$ is the average for each year this
seems fair. 
