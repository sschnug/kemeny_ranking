# Intro
This code solves the [Kemeny Score (Rank Aggregation)](https://en.wikipedia.org/wiki/Kemeny%E2%80%93Young_method) problem exactly. This one
is NP-complete (even for four votes) and therefore hard to do in general.

The approach used is using Integer-programming based on the description of:

> Conitzer, Vincent, Andrew Davenport, and Jayant Kalagnanam. "Improved bounds for computing Kemeny rankings." AAAI. Vol. 6. 2006.

There is also a simple preprocessing based on the *Extended Condorcet criterion* as described in:

> Betzler, Nadja, et al. "Fixed-parameter algorithms for Kemeny scores." Algorithmic Aspects in Information and Management (2008): 60-71.

For some overview of the Integer-programming approach using another solver, see this excellent [blog-post](http://vene.ro/blog/kemeny-young-optimal-rank-aggregation-in-python.html).

This code serves as basic introduction to mathematical-modelling with python, including some non-trivial usage of the scientific-stack (sparse constraint-matrix generation).
It also shows some usage of *cylp*, which is pretty sparse in regards to documentation (and some Cbc-source crawling is needed to achieve good performance here!).

# Requirements / Software used
The code uses:
- python's scientific stack (numpy, scipy)
- [cylp](https://github.com/coin-or/CyLP) as interface to the solver
- [CoinOR Cbc](https://projects.coin-or.org/Cbc) as the (imho best free open-source) core MIP-solver

# Use-cases
See:

> Ali, Alnur, and Marina Meilă. "Experiments with Kemeny ranking: What works when?." Mathematical Social Sciences 64.1 (2012): 28-40.

My use-case: candidate selection in shared appartement castings (although this one is usually feasible to do in a brute-force way).

# Alternatives
- Heuristic: [numerical.recipes C++](http://numerical.recipes/whp/ky/kemenyyoung.html)
- Exact with [FPT-incorporation](https://en.wikipedia.org/wiki/Parameterized_complexity): [kconsens](https://www.akt.tu-berlin.de/menue/software/)
  - Supports commercial and open-source ([GLPK](https://www.gnu.org/software/glpk/)) solvers

# Status
- Prototype-like code (which never produced wrong answers in my tests)
- Not much input-specification / checking
- No support for ties / incomplete votes

# Example
**example_ski_jumping.txt**:

    A: GregorSchlierenzauer SimonAmmann WolfgangLoitzl HarriOlli DimitryVassiliev MartinSchmitt ThomasMorgenstern MartinKoch AndersBardal AndreasKüttel AdamMalysz NoriakiKasai EmmanuelChedal MichaelUhrmann MichaelNeumayer JakubJanda RobertKranjec TomHilde TakanobuOkabe DaikiIto DenisKornilov KamilStoch StephanHocke AndreasKofler JernejDamjan SebastianColloredo StefanHula SigurdPettersen PrimozPeterka PrimozPikl ManuelFettner JanMatura JonAaraas
    B: ThomasMorgenstern GregorSchlierenzauer TomHilde AndersBardal AndreasKüttel SimonAmmann WolfgangLoitzl AdamMalysz AndreasKofler MartinKoch JernejDamjan MichaelNeumayer DimitryVassiliev MartinSchmitt DenisKornilov EmmanuelChedal RobertKranjec HarriOlli MichaelUhrmann ManuelFettner KamilStoch DaikiIto SigurdPettersen NoriakiKasai SebastianColloredo PrimozPeterka JanMatura StephanHocke JakubJanda JonAaraas TakanobuOkabe StefanHula PrimozPikl
    C: AdamMalysz SimonAmmann GregorSchlierenzauer AndreasKüttel ThomasMorgenstern AndreasKofler MichaelUhrmann DimitryVassiliev MartinKoch WolfgangLoitzl AndersBardal MartinSchmitt DenisKornilov TomHilde JakubJanda JernejDamjan SigurdPettersen NoriakiKasai RobertKranjec HarriOlli KamilStoch PrimozPikl ManuelFettner SebastianColloredo MichaelNeumayer TakanobuOkabe JanMatura DaikiIto StefanHula StephanHocke EmmanuelChedal JonAaraas PrimozPeterka
    D: JakubJanda ThomasMorgenstern AndreasKüttel AndreasKofler MichaelUhrmann AdamMalysz WolfgangLoitzl TakanobuOkabe RobertKranjec SimonAmmann MartinKoch MichaelNeumayer JanMatura DimitryVassiliev DaikiIto NoriakiKasai JernejDamjan KamilStoch AndersBardal PrimozPeterka SigurdPettersen ManuelFettner HarriOlli SebastianColloredo MartinSchmitt DenisKornilov EmmanuelChedal StefanHula PrimozPikl StephanHocke TomHilde JonAaraas GregorSchlierenzauer

```python3 run.py example_ski_jumping.txt```

Output:

    Parse input
         ... finished
    Problem statistics
      4 votes
      33 candidates
    Build incidence-matrix
         ... finished
    Solve: build model
      # pairwise constr:  528
        Took 0.000 secs
      # triangle constr:  32736
        Took 0.009 secs
      Extended Condorcet reductions: 230 * 2 relations fixed
    Solve: run MIP

    Welcome to the CBC MILP Solver
    Version: 2.8.12
    Build Date: Feb 22 2016

    command line - ICbcModel -solve -quit (default strategy 1)
    Continuous objective value is 377 - 0.06 seconds
    Cgl0003I 0 fixed, 0 tightened bounds, 1067 strengthened rows, 0 substitutions
    Cgl0004I processed model has 4889 rows, 298 columns (298 integer) and 12446 elements
    Cutoff increment increased from 1e-05 to 1.9999
    Cbc0038I Solution found of 377
    Cbc0038I Before mini branch and bound, 298 integers at bound fixed and 0 continuous
    Cbc0038I Mini branch and bound did not improve solution (0.24 seconds)
    Cbc0038I After 0.24 seconds - Feasibility pump exiting with objective of 377 - took 0.01 seconds
    Cbc0012I Integer solution of 377 found by feasibility pump after 0 iterations and 0 nodes (0.24 seconds)
    Cbc0001I Search completed - best objective 377, took 0 iterations and 0 nodes (0.24 seconds)
    Cbc0035I Maximum depth 0, 0 variables fixed on reduced cost
    Cuts at root node changed objective from 377 to 377
    FlowCoverCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    MIRCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    GomoryCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    ResidualCapacityCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    KnapsackCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    TwoMIRCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    CliqueCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    Probing was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    Gomory was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    Knapsack was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    Clique was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    MixedIntegerRounding2 was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    FlowCover was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)
    TwoMirCuts was tried 0 times and created 0 cuts of which 0 were active after adding rounds of cuts (0.000 seconds)

    Result - Optimal solution found

    Objective value:                377.00000000
    Enumerated nodes:               0
    Total iterations:               0
    Time (CPU seconds):             0.31
    Time (Wallclock seconds):       0.35

    Total time (CPU seconds):       0.31   (Wallclock seconds):       0.35

      CoinOR CBC used 0.350 secs
    Postprocessing
        ... finished
    --------
    SOLUTION
      objective:  377.0
      aggregation:
    ['GregorSchlierenzauer' 'ThomasMorgenstern' 'AndreasKüttel' 'AdamMalysz'
     'SimonAmmann' 'AndreasKofler' 'WolfgangLoitzl' 'DimitryVassiliev'
     'MartinKoch' 'AndersBardal' 'MartinSchmitt' 'MichaelUhrmann'
     'MichaelNeumayer' 'DenisKornilov' 'JakubJanda' 'TomHilde' 'JernejDamjan'
     'RobertKranjec' 'HarriOlli' 'NoriakiKasai' 'EmmanuelChedal'
     'TakanobuOkabe' 'DaikiIto' 'KamilStoch' 'SigurdPettersen' 'ManuelFettner'
     'SebastianColloredo' 'JanMatura' 'PrimozPeterka' 'StefanHula' 'PrimozPikl'
     'StephanHocke' 'JonAaraas']
