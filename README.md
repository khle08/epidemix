# Brief

This package is designed to analyze the network diffusion, such as epidemic spreading and information distribution, motivated by the serious COVID19 pandemic since the end of 2019. The process of network diffusion can be simulated based on the assumption of certain predefined states for each entity and the simulation can be described mathematically under a set of ordinary differential equations (ODE). Except for the models provided by this package, the states and differential equations can also be customized by yourself. In this package, a set of ODE can be applied to the following 2 perspectives:

1. A whole population

2. A Network composed of edges and nodes

Since the virus is the priority task to study, the models provided by this package are all related to epidemics including `SI`, `SIS`, `SIR`, and `SIRV`. In order to analyze the epidemic, the spreading process can be split into the following states so that mathematical model can be built up accordingly:

+ S - Susceptable
+ I - Infected
+ R - Recovered

This is also called SIR model. Sometimes, the model can be even more easier containing only S and I states so that the mathematical mechanism can be better understood. The famous SIR model can not only be applied on the total number of a group of people, but also be able to be implemented to a network composed of nodes and edges. To represent each person by one node and define the relationship between 2 people as the edge connected to the 2 nodes, we can model the society mathematically. A certain epidemic is spreaded throughout the social network from those infected people. By calculating the probabilistic states, i.e:

+ Si(t) probability that node i is `susceptible` at time t 
+ Xi(t) probability that node i is `infected` at time t
+ Ri(t) probability that node i is `recovered` at time t

along with the coefficients such that:
+ β : individual transmission / infection rate
+ γ : recovery rate

We will be capable of tracking the state of each node along the time line. Combining with the undirected graph structure, which is also called a network, the original simple SIR model is transformed from deteministic to probabilistic description.



# Requirements

This package is developed based on the following dependencies:

+ cdlib
+ networkx
+ numpy
+ scipy
+ matplotlib
+ tqdm

p.s. The programming language should be `Python3`.

### Install epidemix

```bash
pip install epidemix
```

And the required dependencies will also be installed automatically.



# The Reproduction Rate

There is a famous factor R0, which is also called: basic reproduction nnumber. This factor indicates the average number of people being infected by an infected person. If R0 > 1, it means that the disease will keep spreading in the society. On the other hand, if R0 < 1, it implies that the infected population will converge and the disease will not spread persistently. 

| Disease    | -            | Transmission     | -        | R0  |
| ---------- | ------------ | ---------------- | -------- | ------ |
| Measles    | 麻疹         | Airborne         | 空气传播 | 12~18 |
| Pertussis  | 百日咳       | Airborne droplet | 空气飞沫 | 12~17  |
| Diptheria  | 白喉         | Saliva           | 唾液     | 6~7    |
| Smallpox   | 天花         | Social Contact   | 社交接触 | 5~7    |
| Polio      | 小儿麻痹     | Fecal-oral route | 粪口     | 5~7    |
| Rubelia    | 风疹         | Airborne droplet | 空气飞沫 | 5~7    |
| Mumps      | 流行性腮腺炎 | Airborne droplet | 空气飞沫 | 4~7    |
| HIV / AIDS | 艾滋病       | Sexual contact   | 性传播   | 2~5    |
| SARS       | 非典型肺炎   | Airborne droplet | 空气飞沫 | 2~5   |



# Simulation Under a Population

The predefined epidemic models can be imported from `epidemix` according to the following code:

```python
from epidemix.macro import SI, SIS, SIR, SIRS, SEIR, SEIRD
```

To activate the model and visualize the results of the model, we need the following extra packages and function called `EpiModel`.

```python
import numpy as np
import matplotlib.pyplot as plt

from epidemix.macro import EpiModel
```

Next, the model can be instantiated by defining the total number of entity in a targeted population, initial infected number, initial recovered number, transmission rate, and recover rate such that:

```python
sir = SIR(1000, I0=50, R0=0, beta=0.4, gamma=0.1)
```

To solve the equation, a period of time should be setup where the interval can be any range according to different conditions.

```python
days = np.linspace(0, 80, 80)
```

Finally, the simulation can be made by `EpiModel`.

```python
epi = EpiModel(sir)
s, i, r = epi.simulate(days)
```

The trend of state transition can be visualized as follows.

```python
fig = plt.figure(facecolor='w')
plt.plot(days, s / sir.N, 'b', alpha=0.5, lw=2, label='Susceptible')
plt.plot(days, i / sir.N, 'r', alpha=0.5, lw=2, label='Infected')
plt.plot(days, r / sir.N, 'g', alpha=0.5, lw=2, label='Recovered')
plt.xlabel('Time /days')
plt.ylabel('Number (1000s)')
plt.grid(True)
plt.legend()
plt.show()
```

![macro.jpg](https://github.com/khle08/epidemix/blob/master/pics/macro.jpg)

## Model Customization - Macro

To customize an epidemic model, we need to write down a set of differential equations in advance.

![sir_macro.jpg](https://github.com/khle08/epidemix/blob/master/pics/sir_macro.jpg)

```python
from epidemix.macro import MacroODE

class SIR(MacroODE):
    def __init__(self, N, I0, R0, beta, gamma):
        self.N = N
        self.I0 = I0
        self.R0 = R0
        self.S0 = N - I0 - R0
        self.initial = (self.S0, I0, R0)

        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, y, t):
        S, I, R = y
        dSdt = -self.beta * S * I / self.N
        dIdt = self.beta * S * I / self.N - self.gamma * I
        dRdt = self.gamma * I
        return dSdt, dIdt, dRdt
```

Mind that the parameters  `self.initial` and `self.N` should be defined in the `__init__` function and the main body of the differential equations should be defined in `derivative` function accordingly.



# Simulation Under a Network

The predefined epidemic models can be imported from `epidemix` according to the following code:

```python
from epidemix.equations import SI, SIS, SIR, SIRV
```

These classes are the default Ordinary Differential Equations (ODE) functions that can be used to simulate in a network. Before starting the simulation, we need the other dependencies, along with the function defined in `epidemix` such that:

```python
import numpy as np
import networkx as nx

from epidemix.epidemic import EpiModel
from utils.plot import draw_probs
```

where `EpiModel` is the most important API being responsible for both network simulation and disease propagation, which shares the same name as the one from `epidemix.macro`. In addition, a given time period is crucial in order to solve ODEs. A timeline should also be generated here.

```python
days = np.arange(0, 10, 0.1)
```

### 1. Network Initialization

Whatever types of network can be generated so that the simulation can be activated based on the network.

```python
num_node = 50
# G = nx.watts_strogatz_graph(num_node, 5, 0.4)     # Small world
# G = nx.powerlaw_cluster_graph(num_node, 5, 0.4)   # Power law
G = nx.gnp_random_graph(num_node, 0.08)             # Random
```

### 2. Instantiation

Take the selected ODEs and Graph (network) into `EpiModel` along with some parameters. Mind that the function will pass `params` into `SIR` ODEs. Namely, the parameters listed here are specifically prepared for SIR function. 

```python
# Note --> SIR  params: I0, R0, beta, gamma
epi = EpiModel(G, SIR, num_state=3, params=[4, 2, 0.4, 0.2])
```

### 3. Simulate

As the parameters are all settled down, the simulation can begin according to the time period. If it is a SIR model, the output would be 3 states where each state is a 2D matrix. The number of row will be defined by the total number of time unit and the number of column will be decided by the total number of node in a network. Each number in the matrix represents the probability that a node staying at THAT corresponding state in a specific moment.

```python
s, i, r = epi.simulate(days)
```

The function will help you get the probability with respect to each time interval. 

![prob.jpg](https://github.com/khle08/epidemix/blob/master/pics/prob.jpg)

### 4. State Propagation

So far, we only get the probabilities of each state for all nodes. However, the deterministic state of each node at time t remains unknown. Although the trend of the probabilities can guide the transformation of each node, we still need to define the sequence first so that the computer can know how to propagate between nodes and between states. In SIR case, S will be turned into I and I will be turned into R.

```pytho
epi.set_propagate(0, 1, neighbor=1, update_state=False)
epi.set_propagate(1, 2, neighbor=None, update_state=False)
status, _ = epi.propagation()
```

`set_propagate` method has 4 parameters (from, to, neighbor, update_state). If SIR is defined properly, 012 will represent SIR respectively and the setup should be done by the number. `neighbor` means that the state transition will happen only when the neighbor of the node has `neighbor` kind of neighbor. S will be infected only when it has $1\rightarrow infected$ neighbors. As for the parameter `update_state`, it is used to deal with the split state transition. If one node can be transformed into 2 optional states, it should follow a sequence. The state that is transformed later should turn it into `True`. 

Finally, the network simulation can be visualized by applying the following function, where status records all the information during network propagation including the actual state of each node, the color of each node, etc. The second parameter indicates what moment we want to observe. The third, forth, and fifth parameters are used to adjust the shape of the plotted result. Therefore, it is better that the number of row and column is in accordance with the number of time interval.

```python
epi.visualize(status, np.arange(16), figsize=(15, 15), n_row=4, n_col=4)
```

![network.jpg](https://github.com/khle08/epidemix/blob/master/pics/network.jpg)



# Self-defined Model: S → I → R

Except for the default epidemic models being defined in `epidemix`, people can also customize their model according to their need. Take SIR model for example here, we assume that the recovered nodes will never get the disease again. The ODE set is formulated as follows:

![sir_eq.jpg](https://github.com/khle08/epidemix/blob/master/pics/sir_eq.jpg)

The adjacent matrix (A) describe the network architecture so that the S nodes can only be contaminated when they have infected neighbors. If there is a connection between 2 nodes, the value would be 1. Otherwise, it would be 0 such that:

```pyth
A = array([[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
           [1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0],
           [0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
           [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0],
           [0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0],
           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1],
           [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1],
           [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0]],
          dtype=int64)
```



# Construct SIR Model with Python Code

The ODE set should be defined in a class inherited from `NetworkODE` class. 

```python
from epidemix.equations import NetworkODE
```

There are also 2 important and identical parts comparing to the one for a pupulation:

1. `__init__` method to initialize the probabilities with respect to different states.
2. `derivative` method to formulate ODE.

```python
class SIR(NetworkODE):
    def __init__(self, A, I0, R0, beta, gamma):
        # numpy 2D Adjacent matrix
        self.A = A
        self.N = len(A)

        # Randomly assign the non-repeated infected and recovered nodes.
        idx = np.random.choice(np.arange(self.N), I0 + R0, replace=False)
        self.I0 = np.zeros((self.N,))
        self.R0 = np.zeros((self.N,))
        self.I0[idx[:I0]] = 1
        self.R0[idx[I0:I0 + R0]] = 1

        # Init matrix should be stacked into a 1D array.
        self.initial = np.hstack([1 - self.I0 - self.R0,    # s(t)
                                  self.I0,                  # i(t)
                                  self.R0])                 # r(t)
        self.beta = beta
        self.gamma = gamma
        self.reproduction_num = beta / gamma    # Definition of "R_0".

    def derivative(self, z, t):
        # The initial "z" starts from "self.initial".
        b = self.beta * z[0:self.N] * np.dot(self.A, z[self.N:2 * self.N])
        r = self.gamma * z[self.N:2 * self.N]
        return np.hstack([-b, b - r, r])
```

If we have 10 nodes in a network, `self.initial` attribute would be a vector with length 10 x \#state​, which is 30 in SIR case. Mind that there are 2 parameters that must be defined here:

1. `self.A` for saving the adjacent matrix.
2. `self.N` for saving the total number of node, which is equal to `len(self.A)`.

As a class is properly defined above, it can be put into `EpiModel` for further simulation. Mind that the parameters defined in the SIR `__init__` class will be set up as `EpiModel` is instantiated with `params` settings.



# Citation

Impact of Vaccination Strategies for Epidemic Node-level SVIR Probabilistic Model. 2020. CL Kuo, MX Chen, WK Victor Chan.



# License

Copyright (c) 2020, Kuo Chun-Lin
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.