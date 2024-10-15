# Passive Imaging

   - [Mitral Cell Mapping](#mitral-cell-mapping)
   - [STIM-ONLY Imaging](#stim-only-imaging)
   - [ODOR-ONLY Imaging](#odor-only-imaging)
   - [ODOR & STIM Imaging - Replicating Mursel’s Results](#odor-stim-imaging)
  
---

## Passive Imaging

_Here, I will give a general overview of the project's work-flow. For specific imaging instructions releted to machinery or software settings, please reference my [IMAGING ROOM](https://docs.google.com/presentation/d/16TTvSaZxpv_2Ob1HF_JzDhX60qj2OGQvm1XAMHGfPwI/edit?usp=sharing) slides._

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

Wavesurfer is used to time the stimulus with respect the the mouse's sniff cycle. This is done as the glomeruli and MCs, in response to identical stimulus parameters, experience sniff-cycle dependent fluctuations i ntheir fluorecence. In this way, while mapping daughter MCs, it is imperative to keep the timing constant of the stimulus with respect to sniff-cycle constant. The latency of the single pulse with respect to sniff-cycle enacted in our protocol is 50 ms. The pulse itself is ___ms long.  

INSERT RUNNING AVERAGE PLOT OF GLOM & MC ACTIVATION WITH RESPECT TO TIMING

For MC mapping, selecting an appropreate power is imperative so that you are not stimulating the surrounding glomeruli beyond what is avoidable. For animals with bright YFP expression at the glomerular layer, apply 2.5V to the AOM which amounts to a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression at the glomerular layer, apply 3.0 to 3.5V to the AOM which amounts to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. 

Use the images produced by MC mapping to further refine the depth at which you would like to image the mitral cells.

**4. Choose your favorite planes and redo MC Mapping**

In the previous step, MC mapping was conducted at an Objective Zoom of 2X. This is so that you would be able to see the activity of a greater area of MCs in response to the glomerular stimulation. However at 2X, the cell bodies of the MCs are not clearly distinguishable. Use the 2X mapping plots to inform the selection of a more precise imaging field where 3 or more MCs are seen to be excited by the stimulation of their respective glomeruli. Further cull this imaging field selection by choosing a panel of ROIs that have minimal cross-talk with neighboring glomeruli. 

Now, you should have a narrowed ROI panel and only a selct few depths that you are interested in. Redo the MC mapping at a 3X zoom while imaging the MCs that you have selected. This step may require you to redraw the ROIs so that they are still located respectively above the MCs which you may have had to move the objective for in 3X. 

_By the end of this section, you should have a select single glomerulus that you would like to stimulate as well as an imaging depth at which you would like to image its repective daughter MCs._

### STIM-ONLY Imaging

In the stimulation-only session, we procede with stimulating the single glomerulus that we chose in the previous step at various latencies. Specifically, the latencies that we have selected are: 10, 30, 60, 120, 180, 240 ms. These latencies were selected for a few reasons. One reason is that these latencies span the entire sniff-cycle, 0 ms being attributed to the mouse's inhilation. In this way, it will be possible to see the fluorescence changes due to the phase of the sniff cycle in response to the same stimulus being presented. Additionally, considering that our aquisition rate is 30 frames/sec, most of these latencies are well divisible by 30, making them easily convertible to frames when analysing the images. 

Stim-only imaging should be conducted in two sessions: one at the glomerul ar level and the next at the MC level. However, it would be best to start with the MC layer as it has been observed that if MCs get stimulated multiple times in the confines of a single day, they diminish significantly in fluorescence. 

The power you should use in stimulating again varies with respect to the ChR2 and GCamp expression in your mouse. For animals with bright expression, apply 2.5V to the AOM which amounts to a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression, apply 3.0 to 3.5V to the AOM which amounts to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. 

Here is an example of some plots that can be rendered from this step. The MATLAB code with which these plots were created is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/Passive_Imaging/SGS_PassiveImg_stimonly_241007.m). 

<img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_ROIs.png" alt="Alt text" width="800"/>

This plot displays the normalized, half-second average fluorescence amplitude of selected regions of interest (ROIs) following stimulation. The different columns correspond to the various stimulation latencies, as indicated on the x-axis. In the top row, which shows the glomeruli within the selected field of view, the stimulated glomerulus at different latencies consistently exhibits the highest average fluorescence amplitude. Similarly, in the mitral cell (MC) layer, the black-outlined dMCs show maximum fluorescence amplitude. 

This plot provides a useful validation of correct glomerulus targeting and shows whether its associated dMCs are co-activated. Since we are using 1P stimulation on a specific glomerulus, the stimulation light has a columnar profile, meaning there's a chance of stimulating axons from other glomeruli passing over the targeted one. Such unintended stimulations can also be identified using this plot.

<img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_dF.png" alt="Alt text" width="800"/>

This plot shows the fluorescence cross-traces of the activation not only of the selected glomerulus but also of its associated dMCs. The fluorescence amplitudes appear to fluctuate in sync with the mouse's breathing cycles, a relationship further clarified by the subsequent average plot.

<img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_stim-time_avg.png" alt="Alt text" width="800"/>

The plot above illustrates the response of the selected glomerulus and its associated dMCs to stimulation at various latencies. As mentioned earlier, despite the constant power of the input stimulation, the responses of both the glomeruli and their corresponding dMCs consistently fluctuate with the mouse’s breathing cycle. Notably, the inhalation phase of the sniff cycle lasts approximately 120 ms, which is reflected in the amplitude of the 180 ms stimulation curve, showing a minimum. Following a full breathing cycle, at 240 ms, the amplitude returns to the maximum observed at the 10 ms stimulation.

This is cool because, unlike odor input, where the volume may vary with the amount of air inhaled, we are directly stimulating the axons of the olfactory sensory neurons (the glomeruli). Yet, even using an artificially generated input stimulus, the response is still modulated by the breathing cycle.

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
