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

<small>Objective Zoom: 1X</small>

Each OB can be modeled using a spere, with glomeruli along the surface and the subsequent excitatory cells extending radially inward. As such, when mapping daughter MCs to their respective parement glomeruli, it is promising to start at a point at the top of this sphere so that the daughter MCs are located directly below your selected imaging plane. This isn't crutial, as at times the OB can have bone overgrowths or simply dim glomeruli, but in any case, this will make your life a tad easier going forth.

Once you think you've found a clear glomerular plane, take a reference image. This will help you to realign the objective if need be. For reference images, it is advised to average 200 imaging frames. 

Next, proceed to explore deeper to see if there are any clusters of MCs that could potentially be imaged. If there are not, go back to the glomerular surface and find a new peice of the bulb to explore.

Repeat this step as many times as necessary.

**2. Take a z-stack under said**

<small>Objective Zoom: 1X</small>

Once you have located your respective plane, take a z-stack or a stack of TIFFs taken at multiple focal distances to create an image with a greater depth of field. Start the z-stack above any glomeruli and set the end of the imaging to below the mitral cells. Take note of your starting and ending depths (Z values) in order to be best able to relocate your favored imaging planes. Choose a step count that would result in increments of 5 µm. MCs have a cell body diameter of 20 to 30 um, but in attempts to find an imaging plane in which the most MCs appear, this increment has proven to be most effective. 

In analysing the resultant z-stack, jot down your favored depths by considering the simple relation below.

$$ \text{favored depth} = \text{starting depth} - (\text{frame number} \times \text{step increment (5µm)}) $$

**3. Conduct MC Mapping at pre-selected depths**

Now that you have a 1X reference image of the glomeruli as well as a list of your favored MC imaging plane depths, it is time to start mapping parent glomeruli to daughter MCs. 

Open the reference image in ImageJ and outline the glomeruli that you would like to stimulate. The light from the stimulation will likely disperse, so make sure to draw the outline within the confines of the glomerulus' perimeter. These selected areas will be refered to as regions of interest or ROIs from now on. Now, save the ROI's that you have selected in a ZIP file and convert them to a format callibrated for the microscope's DMD. 

Using ScanImage, Wavesurfer, and the DMDControl panel, you can now map the MCs. Position the objective (monitoring through ScanImage) to a plane that you would line to image. The stimulation light is columnar and so it is not really affected by the position of the objective along the z-axis. The 2P imaging, however, requires the light to be collected from a focus point and so it is imperative that you position the objective at the imaging plane. 

I would recommend first mapping at the glomerular level (Objective Zoom: 1X) and only then proceed to image at the MC layer (Objective Zoom: 2X) at the selected depths from the prior step. This order will allow you to consider if by stimulating a single glomeruli, you are accidentally also axciting a neighboring one. This can be a result of hitting the axons of a neighboring glomeruli which run along the surface of the OB, provided that the stimulation light is columnar as mentioned previously and doesn't exclusively strike the single intended glomeruli.

Wavesurfer is used to time the stimulus with respect the the mouse's sniff cycle. This is done as the glomeruli and MCs, in response to identical stimulus parameters, experience sniff-cycle dependent fluctuations i ntheir fluorecence. In this way, while mapping daughter MCs, it is imperative to keep the timing constant of the stimulus with respect to sniff-cycle constant. The latency of the single pulse with respect to sniff-cycle enacted in our protocol is 30 µs. The pulse itself is ___µs long.  

[insert running average plot of Glom and MC activation with respect to timing]

For MC mapping, selecting an appropreate power is imperative so that you are not stimulating the surrounding glomeruli beyond what is avoidable. For animals with bright YFP expression at the glomerular layer, apply 2.5V to the AOM which amounts to a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression at the glomerular layer, apply 3.0 to 3.5V to the AOM which amounts to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. 

Use the images produced by MC mapping to further refine the depth at which you would like to image the mitral cells.

**4. Choose your favorite planes and redo MC Mapping**

In the previous step, MC mapping was conducted at an Objective Zoom of 2X. This is so that you would be able to see the activity of a greater area of MCs in response to the glomerular stimulation. However at 2X, the cell bodies of the MCs are not clearly distinguishable. Use the 2X mapping plots to inform the selection of a more precise imaging field where 3 or more MCs are seen to be excited by the stimulation of their respective glomeruli. Further cull this imaging field selection by choosing a panel of ROIs that have minimal cross-talk with neighboring glomeruli. 

Now, you should have a narrowed ROI panel and only a selct few depths that you are interested in. Redo the MC mapping at a 3X zoom while imaging the MCs that you have selected. This step may require you to redraw the ROIs so that they are still located respectively above the MCs which you may have had to move the objective for in 3X. 

_By the end of this section, you should have a select single glomerulus that you would like to stimulate as well as an imaging depth at which you would like to image its repective daughter MCs._

### STIM-ONLY Imaging



### ODOR-ONLY Imaging


### STIM & ODOR Imaging - Replicating Mursel’s Results



---

## Behavioral Training

### Stim-Detection Training

Train the mice to detect the stimulation of their single glom using a Go-NoGo training paradigm.

**Stimulation Power**  
Some mice should start with 3.0V training, others with as high as 3.5V applied to the power source. Mice vary in their channelrhodopsin levels in OSNs. If the detection rate reaches a constant 80%, try increasing the power.

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

**Conditions**  
- AOM set to 3.5V
- DMD has 200 x 200 pixels corresponding to 0.25mm²

**Pulse parameters**  
Training is conducted with 10 pulses. Each pulse is __ ms long with a __ ms break in between pulses.

---

### Lowering Pulse Count

Now, we train the animal to detect shorter pulse sequences until we reach a single pulse.

During a session, it is possible to drop multiple pulse lengths. For example, you can try running behavior for 50 trials with 8 pulses, then 100 at 6, then 150 at 4. Wait until the animal returns to a 90% success rate before decreasing the pulse count.

---

### Psychometric Curve

A psychometric curve will measure the detection rate of the single-pulse stimulation by the power that the stim-laser is set to. The steepest part of this sigmoid curve will be around 75% detection, and we want to understand the power at which this occurs (~80%). 

This helps in two ways:
1. To isolate the principal component of odorant encoding: a single glomerulus and its daughter mitral cells.
2. To measure behavioral output or detection of the stimulation in the presence of a background odorant.

#### Block Trials

Example power-level lists for some block sessions:

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

#### Probe Trial Method (Failed)

We had previously tried making the psychometric curve using non-rewarded probe trials. The probe trials occupied 11% of the Go trials. This small percentage made the acquisition of data very slow considering mice tend to do 300 trials per session. If 150 trials from that session are stimulation Go trials, we would be acquiring approximately 16 data points per session. Furthermore, using probe trials, the mouse would learn the power at which they are mostly rewarded and only lick for the standard power rather than the probes. You would see this on the resultant psych plot day-to-day as the detection became more and more skewed towards the standardized power-level.

INSERT PROBE TRIAL PLOT

## Imaging Session with Behavior

### ODOR & STIM Imaging 

This session should be done on the same day as the ODOR & STIM session that is described below. 

This session is meant to verify minimal cross-glomerular stimulation when we are attempting to target the single glomerulus. 

It is important to do this session prior to the ODOR & STIM session as the odors included in the ODOR & STIM session can linger into the ‘stim only’ trials. If for example you chose an inhibitory odor, the fluorescence changes as a result of stimulation-only can be suppressed by this lingering odorant.

Ensure this session follows after the ODOR & STIM session.

