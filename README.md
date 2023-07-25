# hindawi-paper-data
paper preparation for submission at hindawi limited

The repository consist of:
1) Ns3.36 code with the following extensions:
     1.1) my_dumbbell_sim6.cc in examples/tcp directory - for simulations
     1.2) my_dumbbell module in contrib directory - for creation of dumbbell topology (improvement of Riley code)
     1.3) my_maths_tools module under contrib directory - for mathematical tools including Tikhnov differentiation
     1.4) sim_scripts directory - for simulation scripts as required
     1.5) tcp-mod-newreno.h and tcp-mod-newreno.cc in src/internet/model - for rudementary modification of TcpNewReno as described in the paper
3) Data generated from simulations under generated_data/sim4_tikhonov directory
4) Notebook Python code for analysis and graphs of the results