# Behavioral Training

   - [Stim-Detection Training](#stim-detection-training)
   - [Lowering Pulse Count](#lowering-pulse-count)
   - [Psychometric Curve](#psychometric-curve)
   - [Single-Pulse Stim Detection at Varying Latencies](#single-pulse-stim-detection-at-varying-latencies)
---

_Prior to this training phase, mice were water-deprived until they reached 85% of their original body-weight. This body weight was maintained throughout the following 'Behavioral Training' section._

At the end of this Behavioral Training procedure, the mice should be trained in a Go-NoGo paradigm to lick in response to stimulation of the target glomerulus and to withhold licking when no stimulus is presented. To isolate the principle component of olfactory encoding, we suggest to lower the stimulation to a single pulse. This precision can furhter be instantiated by generating a psychometric curve, which will allow you to determine the threshold power density at which the mouse is no longer able to detect the stimulation.

---
### Stim-Detection Training

The mice were trained to detect optogenetic stimulation of a selected single glomerulus using a Go/No-Go paradigm. Mice were conditioned to lick in response to a presented stimulus and to withhold licking during 'blank' trials where no stimulus was presented. The plot below illustrates the success rates of seven individual mice as they learned the single glomerulus detection task over the course of 6 to 7 days.

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_Training_combinedmice.PNG" alt="Alt text" width="800"/>
</p>

The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_Training_combinedmice.ipynb).

On the x-axis, the sessions are numbered with one session being conducted per day. The length of each session varied in length depending on the motivation of the mouse and the stubborness of the trainer on that day.  Typically, mice began with sessions of 100 to 200 individual trials, gradually working up to sessions with as many as 500 trials. If extending session lengths becomes difficult after a week, the issue may stem from providing too much water per trial.

The y-axis represents the success rate of mice in the detection task. As shown, the detection rate begins at chance level (0.5) on day 1 and steadily increases over time, reaching between 0.9 and 1.0 by the end of training. The shaded region in the plot represents a 95% confidence interval, calculated using the statistical methods outlined [here](https://www.dummies.com/article/academics-the-arts/science/biology/the-confidence-interval-around-a-proportion-149351/).

The most effective stimulus-to-blank trial ratio for both training and continued performance was found to be 0.5 to 0.5. This ratio is optimal because mice tend to have a bias toward licking. If the frequency of stimulus trials exceeds 50%, mice learn to lick on most trials to maximize reward. At a 50% stimulus probability, mice are more motivated to perform the task accurately to avoid the punishment (a 7-second time-out) for licking on a blank trial.

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

Training is conducted with a stimulation consisting of 10 pulses: each pulse having a duration of 10 ms with 6 ms intervals in between.

**Training Using Pavlov**  

**In a Pavlovian conditioning experiment, a contingency is arranged between the presentation of the neutral stimulus, and the delivery of a biologically significant outcome.**  In our case, we stimulate a single glomerulus in the mouse's olfactory bulb and train the mouse to lick a spout for water to indicate detection of the stimulus presented. However, at the onset of training, the naive mouse has neither experienced the stimulus nor learned to lick for water. To fascilitate this association, we employ _pavlovlovian trianing_ to help the mouse connect the stimulus with the reward: _if the mouse detects a stimulus and licks the spout, it will receive water_.

Below is a plot showing the progression of a single mouse's training using Pavlovian conditioning. Initially, the likelihood that a trial will be pavlovian—meaning the mouse receives a reward regardless of its behavior—is quite high. As the training progresses and the mouse learns the contingency between the stimulus and the reward we gradually reduce the frequency of pavlovian trials as seen below.

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_Training_mouse0691.png" alt="Alt text" width="800"/>
</p>

The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_Training_singlemouse.ipynb). 

**Why doesn't the plot showing a combination of all training sessions for all mice exhibit the same check-mark trend?** The single-mouse plot shows a disproportionately higher detection rate during the first few sessions due to the increased likelihood of trials being rewarded via pavlovian conditioning. This doesn’t necessarily indicate that the mouse successfully detected the stimulus, but rather that the reward was automatically given. In the combined training plot, we specifically exclude trials where a Pavlovian reward was administered from the detection-rate calculations. If you review the [code](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_Training_combinedmice.ipynb) for the combined training plot, you'll notice that pavlovian trials were filtered out to provide a more accurate representation of the actual detection performance across all mice.

---
### Lowering Pulse Count

Now we train the animal to detect progressively shorter pulse sequences until it can detect a single pulse.

Within a single training session, multiple pulse lengths can be tested. For instance, you might conduct 50 trials with 8 pulses, 100 trials with 6 pulses, and 150 trials with 4 pulses. While this process is somewhat arbitrary, it is advisable to wait until the animal consistently reaches a 90% success rate before reducing the pulse count further.

Rushing this process can lead to the animal learning to respond to a cue unrelated to the stimulation. In one instance, we reduced the pulse count to a single pulse too quickly and observed a 100% success rate, which seemed unlikely. To investigate, we turned off the stimulation laser, yet the mouse continued to perform at this high detection rate. It became clear that the mouse was responding to an external cue rather than the stimulation itself. We discovered a pressure imbalance at the nose port that only occurred during stimulation trials. After correcting the issue, we had to backtrack and retrain the mouse, which delayed the project by at least two weeks. In summary, it's important to be cautious of unusually high success rates, as they may indicate the animal is learning an unintended cue.

---
### Psychometric Curve

In order to understand at what power to set the stimulation laser to, we need to formulate a psychometric curve. A psychometric curve will measure the detection rate of the single-pulse stimulation as the power that the stim-laser is set to in modulated. Specifically, we would like to understand the stimulation power at which the detection rate is approximately ~80%. Why? The psych curve will be in the shape of a sigmoid function. The y-values of this sigmoid function will range from 50% detection or chance-level to 100% detection. As such, the steepest part of this sigmoid function will be around 75% detection at whatever stim-laser power this detection rate occurs. The steepest part of this curve is the section most sensitive to power shifts and resultant detection changes. As such, it is important to find this detection point for two reasons.

1. We want to isolate the principal component of odorant encoding: a single glomerulus and its daughter mitral cells. In other words, we want to limit the amount of stimulation that hits surrounding glomeruli and would in theory, render detection easier for the animal.
2. In the following section, we will be measuring the behavioral output or detection of the stimulation in the presence of a background odorant. In introducing this odorant, we would like for the stimulus detection to be as sensitive to the changes evoked by the background odor as possible. 


#### Block Method (Succeeded)

After trying various methods of acquiring this psychometric curve, we have concluded on the efficacy of the Block Method. Here, we split a 300 trial session into blocks of 50 trials where 50% of the trials are Go (single-pulse stim) and 50% of the trials are NoGo (no-stim). We pseudo-randomly set the power-level for each block. The first and last blocks should be set to the highest power-level. The power-levels of the blocks in between the first and last block should be randomized so that there aren’t too many blocks at too low of a power (30% of maximum power) where the mouse gives up on detection. An example power-level list for some block-sessions are as follows:

| Day | Maximum AOM Power | Power Levels (% of Maximum Power) |
|-----|---------------------|--------------------------------|
| 1 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 75, 65, 55, 35, 0, 100] |
| 2 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 85, 75, 45, 30, 10, 100] |
| 3 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 95, 60, 50, 40, 15, 100] |
| 4 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 45, 80, 70, 90, 65, 100] |
| 5 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 30, 85, 50, 90, 55, 100] |
| 6 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 80, 65, 95, 20, 60, 100] |
| 7 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 70, 65, 95, 35, 60, 100] |
| 8 | 4.0V applied $58 \ \text{mW}/\text{mm}^2$    | [100, 80, 40, 70, 50, 25, 100] |

Below is an example of a Physichometric Curve taken using the Block Method conducted with only a single mouse. For this specific mouse, only 7 sessions were conducted. This is one of the reasons why the lower power levels have a greater uncertainty. At this point in the project, we only wanted a rough estimate of what the power should be for 80% detection and so this level of data was sufficient.

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_PsychCurve_Block_combinedmice.png" alt="Alt text" width="800"/>
</p>

The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_PsychCurve_Block_combinedmice.py). 

The most important take-away from this method of cunducting the psychometric curve is that the mouse does not learn to lick for only the maximum power-level. This was an issue encountered with the other methods employed as will be discussed below. 

#### Probe Method (Failed)

We also tried rendering the psychometric curve using non-rewarded probe trials. The probe trials were non-rewarded and occupied 11% of the Go trials. We further employed a partial-reward system, where only 80% of the non-probe (maximum stim power) Go trials were rewarded. This was put into place to make sure that the mice didn't learn to lick for only the maximum-power Go trials. Below is a summary figure of the Probe Method for the two mice that we applied the above paradigm with. The shaded region represents a 95% confidence interval which takes into account the number of trials that were aquired for each probe.

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_PsychCurve_Probe_combinedmice.png" alt="Alt text" width="800"/>
</p>

The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_PsychCurve_Probe_combinedmice.py). 

**This data-acquisition method has a few issues which will be discussed below.**

**1. very few data-points collected accross trials**

As mentioned previously, only 11% of the stim trials were probes. This small percentage of probe trials made the acquisition of data very slow. Considering mice tend to do 300 trials per session, if 50% of the trials in that session are stimulation Go trials, we would be acquiring approximately 16 probe data points per session. Below is an example probe-acquisition accross 4 sessions for a single mouse. 

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_PsychCurve_Probe_mouse0070.png" alt="Alt text" width="800"/>
</p>

This plot is mainly used to monitor the data collection process from day to day to see wich power-levels are needed to be probed further for the next session. You can see that for $~60 \ \text{mW}/\text{mm}^2$, the histogram is missing a bar and this is due to 89% (100-11%) of all of the stim trials being conducted at this maximum power density. The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_PsychCurve_Probe_singlemouse.py). 

**2. mouse learns to lick for only maximum power despite partial reinforcement**

Despite imposing a partial-reward paradigm, the mice still managed to learn to lick for the maximal power level at which they were most likely to be rewarded. You can see this on the resultant psych plot as the sessions progressed and the detection progressively became skewed towards the maximum power-level. 

Below is an example plot for the mouse that was analyzed above where the success-rate of the various probes are shown broken into the four seperate sessions run during this data acquisition. 

<p align="center">
  <img src="https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_PsychCurve_Probe_mouse0070_sessions.png" alt="Alt text" width="800"/>
</p>

The code to this plot is linked [HERE](https://github.com/Ekaterina-Koulakova/Single-Glomerular-Stimulation/blob/main/Methods/2_Behavioral_Training/SGS_PsychCurve_Probe_singlemouse_sessions.py). 

In all, not only are you limited in the number of probe-aquisitions per session that you are able to acheive, but you are also limited in the number of sessions that you are able to run before the mouse learns to lick for the standard maximum power.

**3. diffifult to incorporate blank trilas into consideration of stimulation detection rate**

Finally, you may have noticed that on the [Block Method plot](https://github.com/ekaterinakoulak/Single-Glomerular-Stimulation/blob/main/plots/Behavior/SGS_PsychCurve_Block_combinedmice.png), the y-axis limits are 0.5 to 1.0 while for the Probe Method plots above, the y-axis ranges from 0.0 to 1.0. This is because the Block Method allows us to take the coinciding blank trials into consideration. While running probe trials, one can only take a running average of the surrounding blank trials into consideration, but even these blank trials can be successful simply bc the power presented right before the considered range of trials was higher than the one being probed. It is this inability to reliably incorporate blank trials into the detection rate that renders this method weaker. 

Why consider blank trials? A mouse is biased towards licking for water rather than witholding its lick. In other words, a mouse can lick for a stimulation presentation regardless of if it detecting the stimulus or not. For this reason, the mouse witholding its lick on blank trials is just as important to consider and lends creadence to the mouse's decision to eventually lick post-stim presentation.  

#### Randomized Trial Method (Failed)

A plot for this is not unavailable... but it was messy. 

---
### Single-Pulse Stim Detection at Varying Latencies

In this section, we will evaluate whether single-pulse stimulation at the identified 80% detection power can be detected by the mouse at various latencies relative to the sniff cycle. This will be conducted using the Probe Method, where 11% of stimulation trials will be presented at the specified probe latencies below. Stimulations delivered at the default latency will be partially rewarded (80% of trials), while probe trials will not receive any reward.

* Default Latency: 50 ms
* Latency Probes: 10 ms, 30ms, 60 ms, 120 ms, 180 ms, 240ms

**Need to reproduce this plot.**
