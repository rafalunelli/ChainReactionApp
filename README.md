# ChainReactionApp
MATLAB App developed to teach Neutronic Absorption, Nuclear Fission and Chain Reaction. 

Main Concepts: Physical Modelling, Dynamical Systems, ODE Simulations, Collision Theory.

A mathematical model was built in order to transpose the phenomenology of nuclear fission and the chain reaction. The choice of modeling these phenomena is justified by how central these themes are when we try to understand the process of obtaining energy from nuclear reactions. The modeling started from the nuclear reactions that are born from the bombardment with neutrons of the two isotopes of uranium present in the nuclear reactor:

![image](https://user-images.githubusercontent.com/26350626/184990666-1da3b324-60ca-4563-9f9d-1a79522aee4f.png)

The first reaction concerns the nuclear fission of U-235, a process of nuclear transmutation that leads to the formation of lower-mass nuclei and three other neutrons (which eventually can trigger other fissions). At the same time, the second reaction corresponds to the absorption of a neutron by a nucleus of U-238.

An effective way to model interacting and time-varying quantities is through differential equations, which are able to incorporate the parameters and mechanisms of the dynamic system. As we have three entities varying with time, amount of U-235 (U), amount of U-238 (C) and amount of neutrons (n), we get three equations:

![image](https://user-images.githubusercontent.com/26350626/184990685-6d1239f0-dad6-4771-bbea-190570d9a1e0.png)

The equations translate the dynamic behavior of a system whose phenomena are those described by the nuclear reactions mentioned above. The first equation (3) describes the temporal variation of the number of neutrons, which depends directly on nuclear fission (which increases the number of neutrons by two units, if it happens) and on neutronic absorption (by the controlling element, U-238, which absorbs a neutron for each eventual collision with a neutron). It is interesting to note that there is also a proportionality between the magnitude of the temporal variation of neutrons and the other quantities. This is due to the probabilistic character of nuclear reactions. For example, the probability that the reaction will occur depends directly on the amount of elements that can interact. This can also be seen in (4) and (5), where the rate of change of fissile and non-fissile uranium is proportional to the amount of uranium and also of neutrons, that is, of the reacting particles. In these last two equations, the negative sign represents the consumption of the reactants over time. g is a control constant.

![image](https://user-images.githubusercontent.com/26350626/184991575-3e528816-87fb-406f-84dc-d3b339aa486b.png)

(a) Main screen of the developed simulator; (b) Simulation results next to its theoretical prediction by the model.


For questions and suggestion, please send an e-mail to rafaellunelli98 at gmail.com. 

MATLAB Version: R2021a.
