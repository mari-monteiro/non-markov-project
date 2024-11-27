# non-markov-project
This project aims to study the Markovian/Non-Markovian dynamics of a single qubit interacting with a CCA for a homogeneous, disordered and correlated configuration of field modes (local energies). Basically,
we observe a change from non-Markovian to Markovian decay in the presence of disorder by tuning the correlation parameter (\alpha).
 ![ccacerta](https://github.com/user-attachments/assets/362cfaf5-b1cd-4c3f-ab08-66750ac9fa99)
Once the pre print of the article is avalible, I'll upload it in this repository. It may be more efficient then me trying to explain it on this file :). 
https://arxiv.org/abs/2411.14304

Some codes here are quite disturbing, feel free to ask me anything. Maybe one day I'll rewrite it, maybe I won't. I don't know. 
For the fortran code your command line shall look like this: 
marimonteiro@localhost:~$ gfortran -o a.out Markovi5.f90 -L/usr/local/lib -llapack -lblas
Make sure you've installed the lapack lib! 
