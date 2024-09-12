# Table of Contents
1. [Passive Imaging Session](#passive-imaging-session)
   - [Mitral Cell Mapping](#mitral-cell-mapping)
   - [STIM-ONLY Imaging](#stim-only-imaging)
   - [ODOR-ONLY Imaging](#odor-only-imaging)
   - [ODOR & STIM Imaging - Replicating Mursel’s Results](#odor-stim-imaging)
2. [Behavioral Training](#behavioral-training)
   - [Stim-Detection Training](#stim-detection-training)
   - [Lowering Pulse Count](#lowering-pulse-count)
   - [Psychometric Curve](#psychometric-curve)
3. [Imaging Session with Behavior](#imaging-session-with-behavior)
   - [ODOR & STIM Behavioral Imaging](#odor-stim-behavior-imaging)

---

## Passive Imaging Session

_For specific imaging instructions releted to machinery or software settings, please reference my IMAGING ROOM slides. Here. I will only give a general overview of the project's work-flow._

### Mitral Cell Mapping

In order to locate the smallest principle circuit motif responsible for olfactory encoding, we first endevor to locate a glomerulus and its respective daughter mitral cells. 

The mice in which we attempt to do so in are genetically modified to express Channel Rhodopsin in the olfactory sesnsory neurons (OSNs) whose axons terminate in glomeruli on the surface of the olfactory bulb (OB) and GCamp in the respective excitatory projection neurons layers. In this way, we are able to stimulate multiple glomeruli using a DMD pattern light projection and concurrently image at projection neuron layers. Specifically, we image at select depths within the mitral cell (MC) layer, maximizing for cell density. The objective is to locate an imaging plane where at least three daughter MCs can be found and attributed to its stimulated parent glomerulus. At least three cells are necessary in the case that one of the cells is compromised during the subsequent steps of this experiment.

#### Imaging Specifics

**1. Find a clear glomerular plane of interest**

Objective Zoom: 1X

Each OB can be modeled using a spere, with glomeruli along the surface and the subsequent excitatory cells extending radially inward. As such, when mapping daughter MCs to their respective parement glomeruli, it is promising to start at a point at the top of this sphere so that the daughter MCs are located directly below your selected imaging plane. This isn't crutial, as at times the OB can have bone overgrowths or simply dim glomeruli, but in any case, this will make your life a tad easier going forth.

Once you think you've found a clear glomerular plane, take a reference image. This will help you to realign the objective if need be. For reference images, it is advised to average 200 imaging frames. 

Next, proceed to explore deeper to see if there are any clusters of MCs that could potentially be imaged. If there are not, go back to the glomerular surface and find a new peice of the bulb to explore.

Repeat this step as many times as necessary.

**2. Take a z-stack under said**

Objective Zoom: 1X

Once you have located your respective plane, take a z-stack or a stack of TIFFs taken at multiple focal distances to create an image with a greater depth of field. Start the z-stack above any glomeruli and set the end of the imaging to below the mitral cells. Take note of your starting and ending depths (Z values) in order to be best able to relocate your favored imaging planes. Choose a step count that would result in increments of 5 µm. MCs have a cell body diameter of 20 to 30 um, but in attempts to find an imaging plane in which the most MCs appear, this increment has proven to be most effective. 

In analysing the resultant z-stack, jot down your favored depths by considering the simple relation below.

$$ \text{favored depth} = \text{starting depth} - (\text{frame number} \times \text{step increment (5µm)}) $$

**3. Conduct MC Mapping at pre-selected depths**

Now that you have a 1X reference image of the glomeruli as well as a list of your favored MC imaging plane depths, it is time to start mapping parent glomeruli to daughter MCs. 

Open the reference image in ImageJ and outline the glomeruli that you would like to stimulate. The light from the stimulation will likely disperse, so make sure to draw the outline within the confines of the glomerulus' perimeter. These selected areas will be refered to as regions of interest or ROIs from now on. Now, save the ROI's that you have selected in a ZIP file and convert them to a format callibrated for the microscope's DMD. 

Using ScanImage, Wavesurfer, and the DMDControl panel, you can now map the MCs. Position the objective (monitoring through ScanImage) to a plane that you would line to image. The stimulation light is columnar and so it is not really affected by the position of the objective along the z-axis. The 2P imaging, however, requires the light to be collected from a focus point and so it is imperative that you position the objective at the imaging plane. 

I would recommend first mapping at the glomerular level (Objective Zoom: 1X) and only then proceed to image at the MC layer (Objective Zoom: 2X) at the selected depths from the prior step. This order will allow you to consider if by stimulating a single glomeruli, you are accidentally also axciting a neighboring one. This can be a result of hitting the axons of a neighboring glomeruli which run along the surface of the OB, provided that the stimulation light is columnar as mentioned previously and doesn't exclusively strike the single intended glomeruli.

Wavesurfer is used to time the stimulus with respect the the mouse's sniff cycle. This is done as the glomeruli and MCs, in response to identical stimulus parameters, experience sniff-cycle dependent fluctuations i ntheir fluorecence. In this way, while mapping daughter MCs, it is imperative to keep the timing constant of the stimulus with respect to sniff-cycle constant. The latency of the single pulse with respect to sniff-cycle enacted in our protocol is 30 µs. The pulse itself is ___µs long.  

INSERT RUNNING AVERAGE PLOT OF GLOM & MC ACTIVATION WITH RESPECT TO TIMING

For MC mapping, selecting an appropreate power is imperative so that you are not stimulating the surrounding glomeruli beyond what is avoidable. For animals with bright YFP expression at the glomerular layer, apply 2.5V to the AOM which amounts to a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression at the glomerular layer, apply 3.0 to 3.5V to the AOM which amounts to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. 

Use the images produced by MC mapping to further refine the depth at which you would like to image the mitral cells.

**4. Choose your favorite planes and redo MC Mapping**

In the previous step, MC mapping was conducted at an Objective Zoom of 2X. This is so that you would be able to see the activity of a greater area of MCs in response to the glomerular stimulation. However at 2X, the cell bodies of the MCs are not clearly distinguishable. Use the 2X mapping plots to inform the selection of a more precise imaging field where 3 or more MCs are seen to be excited by the stimulation of their respective glomeruli. Further cull this imaging field selection by choosing a panel of ROIs that have minimal cross-talk with neighboring glomeruli. 

Now, you should have a narrowed ROI panel and only a selct few depths that you are interested in. Redo the MC mapping at a 3X zoom while imaging the MCs that you have selected. This step may require you to redraw the ROIs so that they are still located respectively above the MCs which you may have had to move the objective for in 3X. 

_By the end of this section, you should have a select single glomerulus that you would like to stimulate as well as an imaging depth at which you would like to image its repective daughter MCs._

### STIM-ONLY Imaging

In the stimulation-only session, we procede with stimulating the single glomerulus that we chose in the previous step at various latencies. Specifically, the latencies that we have selected are: 10, 30, 60, 120, 180, 240. These latencies were selected for a few reasons. One reason is that these latencies span the entire sniff-cycle, 0 µs being attributed to the mouse's inhilation. In this way, it will be possible to see the fluorescence changes due to the phase of the sniff cycle in response to the same stimulus being presented. Additionally, considering that our aquisition rate is 30 frames/sec, most of these latencies are well divisible by 30, making them easily convertible to frames when analysing the images. 

Stim-only imaging should be conducted in two sessions: one at the glomerul ar level and the next at the MC level. However, it would be best to start with the MC layer as it has been observed that if MCs get stimulated multiple times in the confines of a single day, they diminish significantly in fluorescence. 

The power you should use in stimulating again varies with respect to the ChR2 and GCamp expression in your mouse. For animals with bright expression, apply 2.5V to the AOM which amounts to a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression, apply 3.0 to 3.5V to the AOM which amounts to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. 

Here is an example of some plots that can be rendered from this step.

[INSERT STIM-ONLY PLOTS]

### ODOR-ONLY Imaging

Imaging with an odorant exclusively will allow us to select the one which we will use as our background odorant in subsequent imaging sessions. Here are some potential odorants that are worth screening. 

| Monomolecular Odorant   | Concentration Diluted in <br> 5 µL of Distilled Water |
|-------------------------|-------------------------------------------------------|
| 5Methyl2Hexanone        | 1 µL                                                  |
| Benzaldehyde            | 2 µL                                                  |
| Methyl Valerate         | 5 µL                                                  |
| Ethyl Butyrate          | 8.8 µL                                                |
| Hexanal                 | 10 µL                                                 |
| Heptanal                | 10 µL                                                 |
| 2MBA                    | 10 µL                                                 |
| Propionic Acid          | 20 µL                                                 |

Odor-only imaging will take 4 imaging sessions, imaging not only at the glomeruli and MC planes, but also at 2 different concentrations: 20 and 100 flow. '20 flow' is associated with 1% SVD or saturated vapor pressure and '100 flow' is associated with 5% SVD. 

We are ideally searching for two odorants. 

1. A **strong ligand**

_A **ligand** is an odorant which activates the same olfactory receptors as the single glomeruli which we have selected. In this way, a strong ligand will activate the single glomerulus at even very small concentrations while a weak ligand will only activate the single glomeruli at higher concentrations._

This odorant will activate the single glomeruli in question proceding any of its neighboring glomeruli. In other words, a low-concentration of this odorant should be able to activate this glomeruli. Provided that this project aims to study the circuitry involved in olfactory encoding, this odorant is meant to act as a control for which our selected glomeruli's excitation procedes that of neighboring glomeruli. This in turn will allow the subsequent excitation of our selected glomeruli's daughter MCs to follow [SOURCE].

PROVIDE EXAMPLE PLOT

3. A **weak ligand**

This odorant should activate the selected single glomeruli only after its neighboring glomeruli have been activated. In other words, a low-concentration of this odorant should not be able to activate this glomeruli while a higher concentration should be able to. Provided that this project aims to study the circuitry involved in olfactory encoding, the weak ligand background odorant is meant to provide a source of inhibition for the MCs of our selected glomerulus as it is suspected that the MCs of glomeruli activated prior to our select glomerulus will enact lateral inhibition [SOURCE].

PROVIDE EXAMPLE PLOT

This odorant selection is imperative for the next step 'STIM & ODOR Imaging' where we will explore the interaction between the excitation of the MCs at various latencies with respect to the sniff cycle and background odorant inflicted inhibition.

### STIM & ODOR Imaging - Replicating Mursel’s Results

This session should be done on the same day as the ODOR & STIM session that is described below. 

This session is meant to verify minimal cross-glomerular stimulation when we are attempting to target the single glomerulus. 

It is important to do this session prior to the ODOR & STIM session as the odors included in the ODOR & STIM session can linger into the ‘stim only’ trials. If for example you chose an inhibitory odor, the fluorescence changes as a result of stimulation-only can be suppressed by this lingering odorant.

Ensure this session follows after the ODOR & STIM session.

---

## Behavioral Training

_Prior to this step, it is imperative to water-deprive the mice. We try to keep the mice at 85% of their original body weight prior to water deprivation throughout the steps outlined in this 'Behavioral Training' section._

### Stim-Detection Training

Train the mice to detect the optogenetic stimulation of the single glomerulus that you have selected using a Go-NoGo training paradigm. Mice are expected to lick when there is a stimulus presented and to withhold their lick for 'blank' trials when a stimulus is withheld. The below plot shows the success rates of 7 individual mice learning the single glomerulus detection task over the coarse of 6 to 7 days. 

<img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/SGS_Training_combinedmice.PNG" alt="Alt text" width="800"/>

The code to this plot is linked [here](). 

On the x-axis, you see numbered sessions, one session being conducted per day and varying greatly in length depending on the motivation of the mouse and the stubborness of the trainer on that day. Typically, mice start with sessions that span 100 to 200 individual trials working their way up to sessions spanning up to 500 trials. 

On the y-axis you see the success rate that mice acheive in the detection task. Specifically, you see that this detection rate starts on day 1 at 0.5 or chance level and slowly increases as the days procede to 0.9 to 1.0. The shaded region on this plot shows a 90% confidence interval, calculated using the statistical methods outlined on this [page](https://www.dummies.com/article/academics-the-arts/science/biology/the-confidence-interval-around-a-proportion-149351/).


**Training Stimulation Parameters**  
Mice should start their detection training with power densities between $50$ and $58 \text{mW}/\text{mm}^2$ depending on the strength of their channelrhodopsin (ChR2) expression in the OSNs. If the detection rate curve begins to plateau at a percentage below 90%, it may be necessary to increase your stimulation power. 

| Laser Power Control [V] | Power Density at Window ( $\text{mW}/\text{mm}^2$ )  |
|-------------------------|-----------------------------------|
| 1.47                    | 14                                |
| 1.52                    | 22                                |
| 2.01                    | 22                                |
| 2.50                    | 22                                |
| 3.00                    | 35                                |
| 3.50                    | 50                                |
| 4.00                    | 58                                |
| 4.50                    | 61                                |
_Conditions_
* AOM set to 3.5V
* DMD has 200 x 200 pixels corresponding to 0.25mm²

Training is conducted with a stimulation consisting of 10 pulses: each pulse having a duration of __ ms with __ ms intervals in between.

**Training Using Pavlov**  

In a Pavlovian conditioning experiment, a contingency is arranged between the presentation of the neutral stimulus, and the delivery of a biologically significant outcome, so that the animal learns that a specific stimulus predicts the impending delivery of the unconditioned stimulus.

EXAMPLE PLOT WITH SINGLE MOUSE TRAINING

---

### Lowering Pulse Count

Now, we train the animal to detect shorter pulse-sequences until we reach a single pulse.

During a single training session it is possible to drop multiple pulse-lengths. For example, you can try running behavior for 50 trials with 8 pulses, then 100 at 6, then 150 at 4. This process is a little arbitrary, but I would wait until the animal returns to a 90% success rate before decreasing the pulse count. 

Making the task too difficult too quickly risks the animal learning a cue outside of the stimulation and licking to that cue instead. From experience, we once lowered the pulse count to 1 too quickly and we saw that the mouse was achieving a 100% success rate. This seemed improbable, so we turned off the stimulation laser and the mouse continued to achieve this high detection rate. As such, it became evident that the  mouse was reacting to a cue external to the stimulation itself. It turned out that there was a pressure imbalance at the nose port only during stim trials. Once we fixed this pressure imbalance, we had to revert back and redo some of the training. By trying to speed up the process of dropping the pulse-train length by maybe two days, we cost the project at least an additional 2 weeks of troubleshooting and retraining. More universally, it is beneficial to be suspicious of overly successful behavior as well. 

---

### Psychometric Curve

In order to understand at what power to set the stimulation laser to, we need to formulate a psychometric curve. A psychometric curve will measure the detection rate of the single-pulse stimulation as the power that the stim-laser is set to in modulated. Specifically, we would like to understand the stimulation power at which the detection rate is approximately ~80%. Why? The psych curve will be in the shape of a sigmoid function. The y-values of this sigmoid function will range from 50% detection or chance-level to 100% detection. As such, the steepest part of this sigmoid function will be around 75% detection at whatever stim-laser power this detection rate occurs. The steepest part of this curve is the section most sensitive to power shifts and resultant detection changes. As such, it is important to find this detection point for two reasons.

1. We want to isolate the principal component of odorant encoding: a single glomerulus and its daughter mitral cells. In other words, we want to limit the amount of stimulation that hits surrounding glomeruli and would in theory, render detection easier for the animal.
2. In the following section, we will be measuring the behavioral output or detection of the stimulation in the presence of a background odorant. In introducing this odorant, we would like for the stimulus detection to be as sensitive to the changes evoked by the background odor as possible. 


#### Block Trials

After trying various methods of acquiring this psychometric curve, we have concluded on the efficacy of the Block Trial method. Here, we split a 300 trial session into blocks of 50 trials where 50% of the trials are Go (single-pulse stim) and 50% of the trials are NoGo (no-stim). Here, we manually set the power-level for each block. The 1st and last blocks should be at the highest power-level (4.0V). The power-levels of the blocks in between the 1st and last block should be randomized so that there aren’t too many blocks at too low of a power (30% of maximum power) where the mouse gives up on detection. An example power-level list for some block-sessions are as follows:

| Day | Maximum AOM Power | Power Levels (% of Maximum Power) |
|-----|---------------------|--------------------------------|
| 1 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$   | [100, 75, 65, 55, 35, 0, 100] |
| 2 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 85, 75, 45, 30, 10, 100] |
| 3 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 95, 60, 50, 40, 15, 100] |
| 4 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 45, 80, 70, 90, 65, 100] |
| 5 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 30, 85, 50, 90, 55, 100] |
| 6 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 80, 65, 95, 20, 60, 100] |
| 7 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 70, 65, 95, 35, 60, 100] |
| 8 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 80, 40, 70, 50, 25, 100] |

INSERT BLOCK TRIAL PLOT

#### Randomized Trial Method (Failed)

In this method, we attempted to colloect detection data by purely randomizing the 

INSERT RANDOMIZED TRIAL PLOT

#### Probe Trial Method (Failed)

We also tried rendering the psychometric curve using non-rewarded probe trials. The probe trials were non-rewarded and occupied 11% of the Go trials. We further employed a partial-reward system, where only 80% of the non-probe (maximum stim power) Go trials were rewarded. This was put into place to make sure that the mice didn't learn to lick for only the maximum-power Go trials. This paradigm still experineced a few issues.

For one, this small percentage of probe trials made the acquisition of data very slow. Considering mice tend to do 300 trials per session, if 50% of the trials in that session are stimulation Go trials, we would be acquiring approximately 16 probe data points per session. 

Furthermore, despite imposing a partial-reward paradigm, the mice still managed to learn the maximal power level at which they are most likely to be rewarded. You would see this on the resultant psych plot day-to-day as the detection became more and more skewed towards the standardized power-level.

INSERT PROBE TRIAL PLOT

## Imaging Session with Behavior

Here, we attempt to combine all elements of the experiment up until this point. 

INSERT EXPLANATION OF HYPOTHETICAL PLOTS

INSERT HYPOTHETICAL PLOTS

#### WATER-DEPRIVATION IMPLICATION ON MC IMAGING QUALITY

As you conduct training and take your psychometric curves, it will become evident to you that the quality of your images both at the glomerular and moreso at the MC level become increasingly blurry. Don't panic. This is very likely due to the mouse being in a prolonged state of water deprivation. For this reason, prior to moving onto the next step where we reintroduce imaging, it would be helpful to reintroduce water to the mice. Take some subsequent reference images to verify that water deprivation was in fact the issue with the image quality.

INSERT EXAMPLE WATER DEPRIVATION IMAGES

A few days after reintroducing water, the image quality should significantly improve. Take the water away and shortly after, take a few more days to reinforce the Go-NoGo task with single-pulse stimulation. Once the behavior is back to its original level and the mouse is at 85% of its original body weight, procede to the next imaging step.

### ODOR & STIM Imaging 

We will conduct a Gon-NoGo behavioral session, asking the mouse to lick when it detects a stimulation and to withold its lick when a stimualtion is not present. The trials within the session will be split up between 50% Go trials and 50% NoGo trials. The trial types are as follows:

* stim-only

The first couple (50) of trials in the session will be a stim-only block. This will allow the mouse to not only acclimate to the task at hand, but will also ensure that the stim-only trials are not contaminated by lingering odorants from prior trials.

* empty

Empty trials are those trials for which no stimulation or odorant is issued. These trials will be taken as the NoGo blank trials within the initial stim-only block. This is again to ensure that no lingering odorants contaminate these trials.

* odor-only low concentration

This block will be composed of (20) trials. This will give us 10 trials for which 'odor-only low concentration' images can be taken.

After this block, there will be a block of 20 trials of stim-only which will act as a buffer to reinforce the initial behavioral task.

* odor-only high concentration

This block will be composed of (20) trials. This will give us 10 trials for which 'odor-only low concentration' images can be taken.

After this block, there will be a block of 20 trials of stim-only which will act as a buffer to reinforce the initial behavioral task.

* odor+stim low concentration

This block will be composed of (20) trials. This will give us 10 trials for which 'odor-only low concentration' images can be taken.

After this block, there will be a block of 20 trials of stim-only which will act as a buffer to reinforce the initial behavioral task.
  
* odor+stim high concentration

This block will be composed of (20) trials. This will give us 10 trials for which 'odor-only low concentration' images can be taken.

After this block, there will be a block of 20 trials of stim-only which will act as a buffer to reinforce the initial behavioral task.




