# Behavioral Training

   - [Stim-Detection Training](#stim-detection-training)
   - [Lowering Pulse Count](#lowering-pulse-count)
   - [Psychometric Curve](#psychometric-curve)

---

## Behavioral Training

_Prior to this step, it is imperative to water-deprive the mice. We try to keep the mice at 85% of their original body weight prior to water deprivation throughout the steps outlined in this 'Behavioral Training' section._

### Stim-Detection Training

Train the mice to detect the optogenetic stimulation of the single glomerulus that you have selected using a Go-NoGo training paradigm. Mice are expected to lick when there is a stimulus presented and to withhold their lick for 'blank' trials when a stimulus is withheld. The below plot shows the success rates of 7 individual mice learning the single glomerulus detection task over the coarse of 6 to 7 days. 

<img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/SGS_Training_combinedmice.PNG" alt="Alt text" width="800"/>

The code to this plot is linked [HERE](https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/Methods/Behavioral_Training/SGS_Training_combinedmice.ipynb). 

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


#### Block Trial Method (Succeeded)

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
