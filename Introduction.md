# Introduction

The olfactory system is a relatively simple yet remarkably efficient network dedicated to the detection and processing of a vast array of odor molecules. Despite estimates suggesting there being approximately 30 billion odorous compounds ([Mayhew, 2021](https://www.biorxiv.org/content/10.1101/2020.12.04.412254v2)), the mouse olfactory system encodes these scents using only about 1,000 distinct olfactory receptor types. In doing so, this sensory network filters out irrelevant information to construct a stable percept and decorrelate odor representations.

To study the encoding processes of the olfactory bulb, we focus our attention to the single glomerulus, the principal unit of odor encoding in the olfactory bulb ([Economo et al. 2016](https://doi.org/10.1016/j.neuron.2016.06.001)). Each glomerulus receives input from a specific type of olfactory receptor neuron, and this information is then transmitted to a subsequent population of projection neurons. For now, the projection neuron type which we will be focusing on in this study are mitral cells (MCs). To understand how a signal is transformed from a selected glomeruli to its daughter mitral cells (dMCs), we apply optogenetic stimulation to target individual channelrhodopsin-expressing glomeruli while imaging the activation patterns of subsequent GCamp6f expressing dMCs, concurrently running behavioral assays to link neural responses to mouse perception.

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/projector%20animation.png" width="400">
</p>

Our work building on prior findings which have primarily employed naturalistic odorants as their input stimuli, which result in dynamic activation of the system with poor degrees of glomerular-specificity. For example, work by Dan Rokni showed that as odors of greater similarity to a target odor were added to the mixture presented to the animal, the animal’s detection rate of the target odor dropped at a higher rate than if dissimilar odors were added to the mixture ([Rokni et al., 2014](https://doi.org/10.1038/nn.3775), see Fig. 4). This study suggests that adding odors to the background that have similar glomerular activation profiles as the target degrades the ability of the animal to detect the target odor. This perhaps alludes to a mechanism of lateral inhibition in the region being enacted. Matt Smear et. had similar findings, showing that if a background odorant activates a glomerulus, the detection of optogenetic stimulation targeted at the same glomerulus is substantially diminished ([Smear et al., 2013](https://doi.org/10.1038/nn.3519), see Fig. 2). Both of these studies indirectly lead to the hypothesis that: if neighboring glomeruli are activated, the detection of a local ‘target odor’ or ‘single-glomerular stimulation’ is hindered, alluding to some mechanism of lateral inhibition being enacted.

Building on this prior work, we employ temporally and spatially precise optogenetic stimulation to study the temporal encoding mechanisms within the olfactory bulb. In other words, we study the animals ability to detect a ‘target’ single glomerulus at various latencies with respect to the inhalation cycle.

In 1995, Hopfield first proposed that odor identities are represented through a specific sequence of receptor activations (Hopfield, 1995). The primacy model further refines this concept, suggesting that the initial few glomeruli activated in this encoding sequence have the most significant weight on the odor identity perceived ([Giaffar et al., 2018](https://doi.org/10.1371/journal.pcbi.1012379)). 

<p align="center">
  <img src="https://raw.githubusercontent.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/main/images/mursel_plot.png" width="300">
</p>

Findings from a recent study in our lab, by Mursel Karadas, provides supporting evidence to this primacy model. By stimulating a single glomerulus at varying latencies throughout the sniff cycle, the study observed common excitation of the selected glomerulus’ dMCs throughout the sniff cycle of the animal. Yet, with the introduction of a background odorant, a dynamic response in the dMCs was observed. While the stimulation of the targeted glomerulus consistently evoked an excitatory response at the glomerular level across all stimulation latencies, the related dMCs exhibited an initial excitation early in the sniff cycle followed by inhibition in later phases. This finding would suggest a critical temporal window at the start of the sniff cycle, during which the excitatory signal is able to be transmitted to the cortex. In contrast, later in the sniff cycle, the background odorant appears to suppress this signal, inhibiting cortical transmission.

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/detection_task.jpg" width="300">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/odor_latencies_cropped.png" width="300">
</p>

To test this interpretation of Mursel’s findings while building on the work of Dan Rokni and Matt Smear (aforementioned), I am exploring the implications of learning on this principle olfactory circuit. I am currently running a behavioral detection paradigm, Go-NoGo, where mice are expected to lick when they detect the stimulation and to withhold their lick when the simulation is not presented. Using a DMD, I apply 1P stimulation to excite a single glomerulus at latencies which tile the sniff cycle, alongside a masking background odorant, and simultaneously record the activity of its dMCs using two-photon imaging. 

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/mouse_stim_img.png" width="800">
</p>

There exists several different possible outcomes to this study of single glomerular signal transmission with respect to temporal encoding. Firstly, it is possible that when the odor background suppresses the projection neuron activation at later stages of the sniff cycle, the mouse’s detection will follow. In other words, at the beginning of the sniff cycle as the MC’s are still activated, the mouse will be able to detect the signal, while at later stages when the MC’s activation is suppressed by the odor background, the signal will not be transmitted and the mouse’s detection of the signal will falter. Second, and more interestingly, it may also be possible that even when the MC’s become inhibited at later stages of the sniff cycle, the mouse will be able to detect the stimulation signal. In this case, it warrants to ask, where this signal is being transmitted? If this is the case, we will replicate this experiment while imaging responsive tufted cells (other projection neurons) to see if in the case of behavioral feedback, the signal is transmitted through this secondary type of projection neurons in the olfactory bulb.

Provided that this project is at its inception, I am primarily dealing with designing and testing the experimental procedure to achieve our desired results. In this repository, I will document my progress and any code that I used to analyze my results along the way.

**References**

Mayhew, J, P, O, V. (2021) "Drawing the Borders of Olfactory Space," , undefined(undefined). Available at: https://www.biorxiv.org/content/10.1101/2020.12.04.412254v2.

Economo, N, M., Hansen, K. and Wachowiak, M. (2016) "Control of Mitral/Tufted Cell Output by Selective Inhibition among Olfactory Bulb Glomeruli," Cell Press, 91(2),p. 397-411. Available at: https://doi.org/10.1016/j.neuron.2016.06.001.

Rokni, D. et al. (2014) "An olfactory cocktail party: figure-ground segregation of odorants in rodents," Nature Portfolio, 17(9),p. 1225-1232. Available at: https://doi.org/10.1038/nn.3775.

Smear, C, M. et al. (2013) "Multiple perceptible signals from a single olfactory glomerulus," Nature Portfolio, 16(11),p. 1687-1691. Available at: https://doi.org/10.1038/nn.3519.

Hopfield, J, J. (1995) "Pattern recognition computation using action potential timing for stimulus representation," Nature Portfolio, 376(6535),p. 33-36. Available at: https://doi.org/10.1038/376033a0.

Giaffar, H. et al. (2018) "The primacy model and the structure of olfactory space," Cold Spring Harbor Laboratory, undefined(undefined). Available at: https://doi.org/10.1371/journal.pcbi.1012379.
