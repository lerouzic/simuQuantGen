**Below is the planning for the three days of TD on quantitative genetics:**    

**DAY 1 (December 13th, 2-5 PM)**  

**2-3PM:** Presentation of the simulation tool, brief description of the functions, and some hands-on practice (The Readme.md is a markdown file that will guide you step by step). 
**3-5PM (including 15 minutes break):** Parameters exploration, by binome you will explore the following parameters with the following range (make sure each binome starts the exploration in different order):  
(1) Population size [5,1000];  
(2) Number of loci [1,100];  
(3) Inititial variance [0,4];  
(4) Selection optimum [-10,10];  
(5) Selection strength [Inf, 10-3];
 

**DAY 2 (December 14th, 2-5PM)**  

**2-2:30PM:** Rapid overlook of results compilation and conclusions from Day 1.  
**2:30-3:45PM:** By binome, we will ask you to modify the code of the simulation tool to:  
- (1) incorporate new mutations arising during the simulations, with a mutation rate, mu/individual [10-3;10-2];  
- (2) introduce a selfing rate which may vary, s [0,1];  
- (3) introduce a coefficient of dominance for mutations, h [-2, 2];  
- (4) modify the shape of the fitness function;  
- (5) suppress recombination among loci.  
If you have time you can further explore   
(1)->the correlation;   
(2)->assortative mating (where individuals that are the most alike tend to reproduce amon themselves);  
(3)->different values of h among loci;  
(4)->disruptive fitness function to 1 - exp (-(x-theta)2/2Vs);  
(5)->intermediate values of recombination.  
**3:45-4PM:** Break  
**4-5PM:** Run simulations, start getting graphs  

**DAY 3 (December 15th, 2-5PM)**  

**2-3PM:** Finish simulations and graphs;  
**3-3:45PM:** Interpretations and discussion;  
**3:45-4PM:** Break  
**4-5PM:** prepare slides for friday exam.  
Slides must include a (1) small introduction on the biological question you explored, (2) an explanation on how you modify the code, (3) the results and your interpretation, (4) Conclusion and perspectives.  
