# Passive Imaging

   - [Mitral Cell Mapping](#mitral-cell-mapping)
   - [STIM-ONLY Imaging](#stim-only-imaging)
   - [ODOR-ONLY Imaging](#odor-only-imaging)
   - [ODOR & STIM Imaging](#odor--stim-imaging)
  
---

## Passive Imaging

_This section provides an overview of the workflow for mapping mitral cells. For specific imaging instructions related to equipment or software settings, please refer to my [IMAGING ROOM](https://docs.google.com/presentation/d/16TTvSaZxpv_2Ob1HF_JzDhX60qj2OGQvm1XAMHGfPwI/edit?usp=sharing) slides._

### Mitral Cell Mapping

To identify the principal functional circuit motif involved in olfactory encoding, we begin by locating a glomerulus and its associated daughter mitral cells.

We use genetically modified mice that express Channelrhodopsin in olfactory sensory neurons (OSNs), whose axons terminate in glomeruli on the surface of the olfactory bulb, and GCaMP in the excitatory projection neurons. This allows us to stimulate multiple glomeruli with a DMD light projection while concurrently imaging the projection neuron layers. Specifically, we focus on imaging within the mitral cell (MC) layer at depths that maximize cell density. The goal is to identify an imaging plane where at least three daughter MCs are visible and can be linked to their stimulated parent glomerulus. This number of cells is necessary to account for potential complications in the subsequent steps of the experiment.

#### Imaging Specifics

**1. Identify Field of View**

**Objective Zoom:** 1X

The olfactory bulb (OB) can be modeled as a sphere, with glomeruli located on the surface and the corresponding excitatory cells extending radially inward. When mapping daughter mitral cells to their respective parent glomeruli, it is often useful to start at the top of this sphere, so that the daughter MCs lie directly beneath your selected imaging plane. This isn’t always critical, as the OB can sometimes have bone overgrowths or dim glomeruli, but starting from this position can make your life a tad easier moving forward.

Once you’ve identified a clear glomerular plane, capture a reference image. This will help you realign the objective in future sessions. For reference images, it is recommended to average 200 imaging frames.

Next, explore deeper to check for any clusters of mitral cells that could be imaged. If no clusters are found, return to the glomerular surface and explore a different region of the bulb. 

An example of a well-defined glomerular and associated MC plane is shown below.

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/identify_FoV.jpg" alt="FoV Identification" width="600" />
</p>

This is an interative process so repeat this process as many times as necessary.

**2. Z-stack Imaging**

**Objective Zoom:** 1X

After locating the desired imaging plane, capture a z-stack or a series of TIFF images taken at multiple focal depths to create an image with greater depth of field. Start the z-stack above the glomeruli and set the end of the imaging below the mitral cells. Be sure to record your starting and ending depths (Z values) to facilitate future relocation of your preferred imaging planes. Choose a step increment that results in 5 µm intervals. While mitral cell body diameters range from 20 to 30 µm, a 5 µm increment has been found to be most effective when searching for an imaging plane with maximal MC visibility.

Below is an example z-stack taken at a 2X magnification.

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/240424_MC_Z_stack_2X_15frames.gif?raw=true" width="300"/>
</p>

When analyzing the resulting z-stack, note your preferred depths using the following relation.

$$ \text{favored depth} = \text{starting depth} - (\text{frame number} \times \text{step increment (5µm)}) $$

**3. Conduct MC Mapping**

**Objective Zoom:** 1X (Glom), 2X (MC)

Now that you have a 1X reference image of the glomeruli and a list of preferred mitral cell imaging depths, you can begin mapping the parent glomeruli to their daughter MCs.

Open the reference image in ImageJ and outline the glomeruli you wish to stimulate. Since the stimulation light may disperse, ensure that the outline is drawn within the perimeter of the glomerulus. These selected regions will be referred to as regions of interest (ROIs) from here on. Save the selected ROIs in a ZIP file and callibrate them to match the spacial locations compatible with the microscope’s DMD.

Using ScanImage, Wavesurfer, and the DMDControl panel, you can proceed with MC mapping. Position the objective (viewing through ScanImage) to the imaging plane you wish to focus on. The stimulation light is columnar, so it is not significantly affected by changes along the z-axis. However, 2P imaging requires a focused light collection, so it's critical to position the objective at the desired imaging plane.

I recommend first performing the mapping at the glomerular level (Objective Zoom: 1X), and then imaging at the MC layer (Objective Zoom: 2X) at the selected depths from the previous step. This order allows you to assess whether stimulating a single glomerulus inadvertently excites a neighboring one. This could occur if the stimulation light, being columnar, stimulates axons from nearby glomeruli that run along the surface of the olfactory bulb (OB).

<!--
Wavesurfer is used to synchronize the stimulus with the mouse's sniff cycle, as both glomeruli and MCs experience sniff-cycle-dependent fluorescence fluctuations in response to identical stimuli. It’s crucial to maintain constant stimulus timing relative to the sniff cycle while mapping daughter MCs. The protocol sets a latency of 50 ms between the stimulus pulse and the sniff cycle, with a pulse duration of 10 ms.
-->

For MC mapping, it’s important to choose the correct stimulus power to minimize the stimulation of surrounding glomeruli. For animals with bright YFP expression in the glomerular layer, apply 2.5V to the AOM, which corresponds to a $22 \ \text{mW}/\text{mm}^2$ power density at the objective surface. For animals with dim YFP expression, apply 3.0 to 3.5V to the AOM, resulting in a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density at the objective surface.

Use the images obtained from MC mapping to further refine the depth at which you wish to image the mitral cells.

**4. Refine Imaging Planes**

**Objective Zoom:** 1X (Glom), 3X (MC)

In the previous step, MC mapping was conducted at an Objective Zoom of 2X to capture a larger area of MC activity in response to glomerular stimulation. However, at 2X, the cell bodies of the MCs are not clearly distinguishable. Use the 2X mapping results to identify a more precise imaging field, where at least three MCs are excited by stimulation of their respective glomeruli. Further refine this field by selecting ROIs that minimize cross-talk with neighboring glomeruli.

Now that you have a narrowed set of ROIs and selected depths of interest, conduct a second round of MC mapping at a 3X zoom to focus on the MCs you’ve identified. This step may require you to redraw the ROIs to ensure they are still positioned correctly above the MCs, especially if you’ve had to adjust the objective to accommodate the 3X zoom.

Below are example images from a successful MC mapping session at depths of 0 µm and 92 µm. The red region in the left-most column of images indicates where the stimulation was applied. The right-most bottom image shows the four daughter MCs that responded to the single-glomerular stimulation.

<p align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/images/MC_Mapping.jpg" alt="FoV Identification" width="800" />
</p>

_By the end of this section, you should have selected a single glomerulus to stimulate and identified the imaging depth at which to capture its respective daughter MCs._

### STIM-ONLY Imaging

In the previous section, we mapped the daughter MCs to the glomerulus without synchronizing the stimulus to the sniff cycle. As a result, stimulation was presented independently of the sniff cycle. However, the activation of daughter MCs is highly dependent on the phase of the sniff cycle during which they are stimulated.

In the STIM-ONLY session, we proceed by stimulating the single glomerulus at various latencies relative to the sniff cycle. The selected latencies are: 10, 30, 60, 120, 180, and 240 ms. These latencies were chosen for a few reasons. Firstly, they span the entire sniff cycle, with 0 ms corresponding to the mouse's inhalation phase. This allows us to observe fluorescence changes in response to the same stimulus presented at different phases of the sniff cycle. Additionally, given our acquisition rate of 30 frames per second, these latencies are well-divisible by 30, making it easy to convert them into frame numbers when analyzing the images.

STIM-ONLY imaging should be conducted in two sessions: one at the glomerular level and another at the MC level. However, it is recommended to start with the MC layer, as repeated stimulation of MCs within the same session can significantly diminish their fluorescence over time.

The power used for stimulation varies depending on the expression levels of ChR2 and GCaMP in the mouse. For animals with bright expression, apply 2.5V to the AOM, which results in a $22 \ \text{mW}/\text{mm}^2$ power density at the surface of the objective. For animals with dim YFP expression, apply 3.0 to 3.5V to the AOM, corresponding to a $35$ to $50 \ \text{mW}/\text{mm}^2$ power density.

Example plots generated from this step can be found using this [MATLAB code](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/1_Passive_Imaging/SGS_PassiveImg_stimonly_241029.m). 

<div align="center">
   <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_ROIs.png" alt="Alt text" width="2000"/>
</div>

This plot displays the normalized, half-second average fluorescence amplitude of selected ROIs following stimulation. The fluorescence delta from the no-stim control condition is what is shown as this condition was subtracted from each image. The different columns correspond to the various stimulation latencies, as indicated on the x-axis. In the top row, which shows the glomeruli within the selected field of view, the stimulated glomerulus at different latencies consistently exhibits the highest fluorescence amplitude. Similarly, in the MC layer, the black-outlined daughter MCs show maximum fluorescence amplitude.

This plot helps validate correct glomerulus targeting and checks whether its associated daughter MCs are co-activated. Since we use 1P stimulation on a specific glomerulus, the stimulation light has a columnar profile, meaning that axons from neighboring glomeruli may also be stimulated. Such unintended stimulations can be detected through this plot.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_dF.png" alt="Alt text" width="800"/>
</div>

This plot shows the fluorescence cross-traces of the selected glomerulus and its associated daughter MCs. The MC fluorescence amplitudes fluctuate in sync with the mouse's breathing cycles, a relationship further clarified by the subsequent averaging plot.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_stim-time_avg.png" alt="Alt text" width="800"/>
</div>

The plot above illustrates the half-second-averaged response of the selected glomerulus and its associated daughter MCs integrated after stimulation presentation at various latencies. Despite constant power in the input stimulus, the responses of both the glomeruli and their corresponding MCs fluctuate with the mouse’s breathing cycle. Notably, the inhalation phase of the sniff cycle lasts approximately 120 ms, which is reflected in the amplitude of the 180 ms stimulation curve, showing a minimum. After a full breathing cycle, at 240 ms, the amplitude returns to its maximum value, similar to the response at 10 ms.

This is cool because, unlike odor input, where the volume may vary with the amount of air inhaled, we are directly stimulating the axons of the olfactory sensory neurons (the glomeruli). Yet, even with an artificially generated input stimulus, the response is still modulated by the breathing cycle.

### ODOR-ONLY Imaging

Imaging with a passively-presented odorant will allow us to select a background odorant for subsequent imaging sessions. Below are some potential odorants worth screening which activate many dorsal glomeruli.

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

Odor-only imaging will require four imaging sessions, covering both the glomeruli and MC planes, at two different concentrations: 20 and 100 flow. The "20 flow" corresponds to 1% SVD (saturated vapor density), while the "100 flow" corresponds to 5% SVD.

We are ideally searching for two odorants described below. 

1. A **strong ligand**

_A **ligand** is an odorant which activates the same olfactory receptors as the single glomeruli which we have selected. In this way, a strong ligand will activate the single glomerulus at even very small concentrations while a weak ligand will only activate the single glomeruli at higher concentrations._

A **strong ligand** will activate the single glomeruli in question proceding any of its neighboring glomeruli. In other words, a low-concentration of this odorant should be able to activate this glomeruli. Provided that this project focuses on olfactory encoding circuitry, this odorant will serve as a control, ensuring that the excitation of the selected glomerulus precedes the activation of adjacent glomeruli, allowing for the subsequent activation of the daughter MCs.

3. A **weak ligand**

A **weak ligand** will activate the selected single glomeruli only after its neighboring glomeruli have been activated. In other words, a low-concentration of this odorant should not be able to activate our selected glomerulus while a higher concentration should be able to. Provided that this project aims to study the circuitry involved in olfactory encoding, the weak ligand background odorant is meant to provide a source of inhibition for the MCs of our selected glomerulus as it is suspected that the glomeruli activated prior to our select glomerulus will enact lateral inhibition.

This odorant selection is crucial for the next step, STIM & ODOR Imaging, where we will explore the interaction between MC excitation at various latencies relative to the sniff cycle and the inhibition caused by the background odorant selected.

Here are example plots generated during the odor screening process using this [MATLAB code](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/1_Passive_Imaging/SGS_PassiveImg_odoronly_241029.m). 

The plot below shows the normalized, half-second average fluorescence amplitude of selected ROIs following inhalation onset, one second after odor (Hexanal) presentation. The fluorescence delta from the 'no-odor' control condition is shown, as this condition was subtracted from each image. All ROI activities were normalized relative to their respective anatomical layers before averaging. The columns represent two different odor concentrations (1% and 5% SVD), indicated on the x-axis. In the top row, which displays the glomeruli within the selected field of view, it can be seen that the selected glomerulus is barely affected by the odor at the lower concentration, while at the higher concentration, the glomerulus shows slight activation. The black-outlined daughter mitral cells (MCs) exhibit the same trend.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/odor_ODOR-ONLY_ROIs_Hexanal.png" alt="Alt text" width="800"/>
</div>

Below are the fluorescent cross-traces of the selected ROIs, normalized with respect to their anatomical layer.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/odor_ODOR-ONLY_dF_Hexanal.png" alt="Alt text" width="800"/>
</div>

In the glomerular layer, I plotted the averaged cross-traces of the selected glomerulus’ activation at two concentrations (in blue). The black line represents the 'no-odor' condition. As shown, the fluorescence cross-trace of the 'no-odor' condition mirrors that of the odor conditions, indicating cross-contamination between sessions. This suggests that odor lingers in the nose port and rubber tubes used for odor delivery between trials. To avoid this issue, it may be helpful to conduct passive odor imaging using the following block structure:

- **1st block:** 10 trials at the no-odor condition
- **2nd block:** 10 trails at a low-concentration odor condition
- **3rd block:** 10 trails at a high-concentration odor condition

In the glomerular plot, cross-traces from neighboring glomeruli activated by odor presentation are also visible. The neighboring glomeruli’s activation precedes that of our selected glomerulus, which is indicative of a weak-ligand odor.  However, ideally, daughter MCs would be inhibited in the case of a weak-ligand odor, suggesting lateral connectivity between neighboring and selected glomeruli. Here, the activity of the daughter MCs correlate with that of the parent glomerulus, indicating low lateral connectivity.

Although no stimulation is presented in the **ODOR-ONLY** condition, we are working towards combining stimulus and odor presentations. The following plot was created to mirror the time-averaged response in the **STIM-ONLY** condition. 

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/odor_ODOR-ONLY_stim-time_avg_Hexanal.png" alt="Alt text" width="800"/>
</div>

Here, you can see the response of the selected glomerulus and its associated daughter MCs at various latencies before stimulation is presented. The plotted dots represent the half-second time average of the fluorescence trace after each time point. 

Overall, **these **ODOR-ONLY** plots are not examples of an ideal odor for this experiment**, but rather a starting point for odor screening. The ideal odor would be a weak ligand that activates the targeted glomerulus only at a higher concentration, after neighboring glomeruli have been activated. Additionally, to ensure connectivity between neighboring glomeruli and the selected glomerulus, it is important that the activity of the daughter MCs does not perfectly correlate with the activity of the selected glomerulus. Ideally, the signal should be modified in some way — for example, even though the selected glomerulus is excited, the daughter MCs might be inhibited.

## ODOR & STIM Imaging 

Here are two of the weaknesses carrying over from the previous **ODOR_ONLY** procedure in odor selection (Hexanal). 
1. The odor: _Hexanal_, that we used did not demonstrate appropreate connectivity with its neighboring glomeruli.
2. We witnessed contamination between trials, meaning the odor would linger in the nose port for longer than 15 sec.

Given these weaknesses in the odor which we selected to procede with, the subsequent plots that were generated to analyze the selected glomerulus' activation in the **ODOR-STIM** condition did not turn out as expected. Below I will explain more in depth the issues with this iteration of the experimenet. 

Alas, this is the code that I wrote to analyze the results of this session: [MATLAB code](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/1_Passive_Imaging/SGS_PassiveImg_stimodor_241102.m).

The **ODOR-STIM** session is very long. Below I summarized the trials presented for each condition in a randomized order.

- **no stim nor odor:** 10 trials
- **odor-only high concentration:** 10 trials
- **odor-only low concentration:** 10 trials
- **stim-only:** 10 trials per 6 latencies
- **stim & odor low concentration:** 10 trials per 6 latencies
- **stim & odor high concentration:** 10 trials per 6 latencies

**DEDUCTION FROM LAST ITERATION OF PROJECT**

The fact that their is cross-contamination of the odor from trial to trail was present even in this paradigm. To avoid this, in the next iteration, we will break this session into 3 seperate sessions.

- **1st session :** stim-only and no-odor trials
- **2nd session :** stim-only, no-odor, odor-only low concentration, stim & odor low concentration
- **3rd session :** stim-only, no-odor, odor-only high concentration, stim & odor high concentration

By breaking the one massive session into three seperate sessions and presenting different condition as blocks of trials rather than randomized, we will not only avoid contamination in the no-odor and stim-only sessions, we will have the no-odor and stim-only trials serve as a control in subsequent sessions where odor is introduced.

The below plot shows the cross-traces of the **STIM-ONLY** condition at the 6 various latencies. You can see that relative to the plots generated in the **STIM-ONLY** section, [LINK to STIM-ONLY plot](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_dF.png), these traces appear to be contaminated with an odor as the activation of the selected glomeruli is relatively suppressed. 

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_STIM-ONLY_dF.png" alt="Alt text" width="800"/>
</div>

Below is the half-second stim-time averaged plot from the **STIM-ONLY** condition. You can see more clearly here that despite the activity of the glomeruli being excitatory, the activity of the MCs are suppressed. This clearly demonstrates the presence of odor contamination. Please reference the previous **STIM-ONLY** plot to compare these results: LINK to STIM-ONLY plot](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stim_STIM-ONLY_stim-time_avg.png).

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_STIM-ONLY_stim-time_avg.png" alt="Alt text" width="800"/>
</div>

Below are fluorescent cross-traces from the glomerular and mitral cell layers under an **ODOR-ONLY** condition. These cross-traces are successful.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_ODOR-ONLY_dF_Hexanal.png" alt="Alt text" width="800"/>
</div>

Below are half-second stim-time average curves from both the **STIM-ONLY** and **ODOR-ONLY** conditions from the larger session.

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_ODOR-ONLY_STIM-ONLY_stim-time_avg_Hexanal.png" alt="Alt text" width="800"/>
</div>

Now we address the **STIM-ODOR** condition. Below are the fluorescence cross-traces from the two odor concentrations (columns) for each cellular-layer (rows). You can see that the stimulation administed at various latencies barely affects the odor & stim-curves. As such, it can be concluded htat the odor concentrations which we are using in this step are too high. This introduces another amendment to the next iteration of the project.

**Prior to imaging the ODOR & STIM condition** we need to conduct a psychometric curve on the minimal power necessary for the more to detect the stimulus and the minimal odor-concentration necessary to mask the detection of the stimuls. 

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_STIM-ODOR_dF_Hexanal.png" alt="Alt text" width="800"/>
</div>

Finally, below is a compilation of all of the half-second stim-time averages from each condition. This si a confusing plot and is in the process of being simplified. However, it is a good analysis for the experimenter to use to see what is going on. 

<div align="center">
  <img src="https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/plots/Imaging/stimodor_STIM-ODOR_stim-time_avg_Hexanal.png" alt="Alt text" width="800"/>
</div>

Overall, despite the plot above being generated with an odor which is not conducive for this project, they are shown to demonstrate the availability of code to analyze results after **STIM-ONLY**, **ODOR-ONLY**, and **STIM & ODOR** sessions. 

**SUMMARY of ODOR & STIM SESSION TAKE-AWAYS FOR THE NEXT SESSION**

- select a weak-ligand odor which has high connectivity with neighboring glomeruli which procede in activation
- break the large conglomerate session into smaller sessions to avoid contamination
- redo this passive imaging session after behavioral training to select odor-concentration which are conducive to the stimulation power being used ie. the high-concentration of the odor should be able to 'mask' the stimulation.
