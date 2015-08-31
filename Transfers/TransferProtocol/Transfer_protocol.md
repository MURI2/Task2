## MURI2 Task 2 Sampling Protocol

**Note:** *I am still working out small details (getting undergrad help for prep) so this is incomplete.*

**Note:** *Work is split into 10-day units. This is so we can guarantee that we have enough supplies for 10-day transfers.*

### Goal: To ensure transparency in Task2

### Shaker layout

There are two Innova 44 shakers in the Lennon Lab. Each has a maximum capacity of 96 50 mL Erlenmeyer flasks.

The top shaker lines that are transferred at a high volume (selection) and is marked with the following: 

<span style="font-size:1em;">\\[\left | N_{e}*s \right |> 1.0\\]
</span>

The bottom shaker contains lines that are transferred at a low volume (drift) and is marked with the following:

<span style="font-size:1em;">\\[\left | N_{e}*s \right |< 1.0\\]
</span>

Within each shaker the lines are grouped into five blocks. 

To briefly go over the layout, we have two transfer sizes, three transfer times, 6 strains, and replicates of each treatment. 2 * 3 * 6 * 5 = 180. So within each shaker there are 90 lines, all of those lines are transferred at either the high or low dilution are contained within a shaker. The 90 lines are divided into 5 blocks (the number of replicates) and within each block the position of transfer time and strain permutations are randomized using a custom Python script that is on GitHub.

The layout of the top shaker is as follows:
![](../Images/TopShakerFormat.png =500x300) 

The layout of the bottom shaker is as follows:
![](../Images/BottomShakerFormat.png =500x300) 

I will be adding borders indicating blocks for each 96-flask shaker tray soon and adding stickers to the bottom of each flask clamp to indicate position on a 96-well format.

### Labelling system

Our labelling system differs from that used in Task 1, primarily because of the different experimental design. The primary goal in creating this label system was to place all relevant metadata in a label that is the same length across all variations of the label. 

The label consists of 6 elements, 3 of which being numerical and 3 being alpha-numerical. 

A table describing the order and element contents is below:

| Position    | Describes | Elements |
|:------------:|:--------------------------------:|:-------------------------:|
| 1       |                Transfer time               |           0, 1, or 2          | 
| 2        |               Shaker              |           T (top) or B(bottom)          |           
| 3    |              96-well position             |           A01 - H12           |          
| 4     |              Dilution              |           D1, D3, or D5          |              
| 5 |              Strain             |            A, B, C, D, F, or J           |             
| 6 |              Replicate (block)              |           1-5         | 
 

Below are tables describing what each element means in the labelling system that is not intuitively clear.

##### Element 1: Transfer time.

| Transfer time | Description |
|:------------:|:--------------------------------:|
|0 | 2-Day Transfers|
|1 | 10-Day Transfers|
|2 | 100-Day Transfers|

##### Element 4: Dilution.

| Abbreviation | Description |
|:------------:|:--------------------------------:|
|D1 | 1:10 dilution|
|D3 | 1:1,000 dilution|
|D5 | 1:100,000 dilution|


##### Element 5: Strain
| Abbreviation | Description |
|:------------:|:--------------------------------:|
|A | *Arthrobacter sp.* KBS0703
|B | *Bacillus subtilis* 168
|C | *Caulobacter crescentus* NA1000
|D | *Deinococcus radiodurans* BAA-816
|F | *Flavobacterium sp.* KBS0721
|J | *Janthinobacterium sp.* KBS0711

##### Element 6: Replicates
| Abbreviation | Description |
|:------------:|:--------------------------------:|
|1| Block 1 (top-left)
|2| Block 2 (top-right)
|3| Block 3 (bottom-left)
|4| Block 4 (bottom-right)
|5| Block 5 (center)


#### Transfer time

The color of the sticker which the label is printed on indicates the transfer time of that line.


2-Day transfers            |  10-Day transfers         | 100-Day transfers
:-------------------------:|:-------------------------:|:-------------------------:|
![](../Images/greendot.png =150x150)  |  ![](../Images/yellowdot2.png =140x140)  |  ![](../Images/reddot.jpg =150x150)

##### Example

Say we have a line of Janthinobacterium that is being transferred every two days at a high transfer volume. It is the fifth replicate and it has been randomly assigned to position G10 on the 96-well format shaker tray. 

Our label would be **0-T-G10-D1-J-5**. 

### Prep

Every 10 days we will need a minimum of 360 50 mL Erlenmeyer flasks filled with 10mL of media. During 10-day increments where there is a 100-day transfer we need 420 flasks. It is good practice to make sure that there are at least ~ 20% more flasks with media at hand.

#### Media

The recipe for the complex media created for Task 2 is PYE with 0.2% glucose and 0.1% casamino acids. The recipe is as follows: 

+ 2 g bactopeptone
+ 1 g yeast extract
+ 0.3 g MgSO<sub>4</sub> x 7 H<sub>2</sub>O
+ 2 g dextrose
+ 1 g casamino acids

Bring up to 1 L final volume with DI H<sub>2</sub>O
Add 15 g Agar/ L if you want plates
Sterilize by autoclaving

#### Flasks

Labelled flasks for 10 days of transfers (5 2-day transfers and 1 10-day transfer) will be prepped every 10 days. 

For prepping flasks, the protocol is as follows. Using a calibrated liquid dispenser, aliquot 10 mL of the MURI media into clean 50 mL Erlenmeyer flasks. Place a 30 mL plastic beaker over each flask. Then, sterilize a tray of those flasks marked with autoclave tape with ~3/4 in. of water for 30 minutes in the autoclave. Once you've removed the flasks and they are cool, place them onto a dry tray and dump the water. 

#### Tubes

The same outline for flasks applies to tubes. Fill with 10 mL of media, autoclave with metal caps. 

#### Tips

We are using USA Scientific filtered pipet tips for our transfers. The Lennon Lab has a set of hood-dedicated pipettes, so those will always be available for transfers. Will will be ordering and keeping track of pipet tips over the course of the experiment. 

#### Printing labels

The IP address of the printer in the Lennon Lab is on the printer. To print the stickers, do the following:

+ Choose sticker color by what transfer times you're printing (see above image). These stickers will be above the printer in the Lennon Lab.

+ Open up the sticker format in the Box folder ~/MURI2/Aim2/Transfers/Labels. To get all 60 labels for a particular transfer time, you'll need to print three sheets of stickers.

+ Open the slot on the printer marked with orange tape and insert your sticker sheet sticker side up.

+ Print the stickers using the word document. 

For a normal two day transfer you will need 60 stickers. 

#### Setup 

Ten days of flasks will be labelled and organized before the start of a ten-day cycle. In each ten day cycle there are five two-day transfers and one ten-day transfer, so there a total of 6 trays that need to be prepped in a normal week (not including weeks where there are 100 day transfers). Each tray is labelled with a sticker at the position where the flask should be placed.  

#### Washing 

Dump the contents of each used 50 mL flask or tube into the large flask next to the sink in the Lennon Lab. 

Fill each flask half-way or tube with water and autoclave to help remove any residual biofilms, then proceed with washing.

For washing the 360 plus 50 mL flasks that will be dirtied each week, we've purchased a 81 IXC 14 dishwasher rack from Lancer. This rack is compatible with the fourth floor dishwasher and is capable 

We are waiting on the rack. Once it arrives we will find what protocol is best for the fourth floor dishwasher. 10 mL test tubes can be washed in the tray that's in the dishwasher (~140 tubes/ 10 days). 

### Transfers

I have created a Google calendar detailing when transfers are occurring and it is available [here](https://www.google.com/calendar/embed?src=7cvon84l098k62h7aqt7gdepqg%40group.calendar.google.com&ctz=America/New_York )

#### Transfer procedure

When someone walks in to do a normal two-day transfer there should be a tray already prepared. Each tray is divided into two sections. One for the top shaker and one for the bottom shaker. Within each section there are five rows of six. Each row corresponds to a block on the 96-well formatted tray and tells you what you are transferring that day. Double check to make sure that the flask order is correct.

Place your pipet tips and the vortexer in the hood and UV for ~10-15 minutes. Press the bottom on the hood that powers the outlet inside the hood.

Starting with the top shaker, remove the tray and then turn on the UV light. Place the lines you're transferring that day on a new tray that has been wiped down with ethanol. Bring your tubes to the hood and ethanol your gloves. 

Working across each block, transfer 1 mL of culture from the old flasks into the new ones.

Carry the tray of newly transferred lines over to the 96-well formatted tray and place each line into its 96-well position, indicated by the third object on the label on the sticker. Turn off the UV light and carefully place the tray back into the shaker.

Do the same thing for the bottom shaker with the exception of dilutions. 

For lines with "D3", you only need to transfer 10 uL into the new flasks (a 1:1,000 dilution). For lines marked with "D5" you need to transfer 10 uL into a 10 mL tube, briefly vortex at high speed, and aliquot 100 uL into the 50 mL flask (a 1:100,000 dilution). 

Take all the old flasks and put them in a tray labelled with today's date. Take the flasks from the previous transfer and label them as waste. They'll be cleaned by (insert undergrad's name here).

#### Human error

A goal with the protocol is to minimize the probability of mistakes occurring. However, if at any point you doubt the sterility of some supply during a transfer (ex. pipet tip brushes against glove, 10 mL dilution tube splashes out, etc), get a new one. There are always spares ready at hand.

If you notice that a flask is translucent (i.e. no growth), go back to the lines from either two or four days ago and transfer from that flask. Record the date, what line didn't grow, and how many steps back you had to transfer from (1 for two-days back, 2 for four-days back if the issue was with a two-day transfer). 

### Responsibilities

Will Shoemaker 

+ Overseeing Task2
+ Transfers (~1/week)
+ Ordering
+ Collecting samples
+ generating and analyzing data

Emily Williams

+ Transfers (~1/week)
+ Collecting samples

Sam Miller
 
+ Transfers (~1/week) 
+ Collecting samples

yet-to-be-named-undergrad 

+ Making media 
+ Prepping flasks
+ Waste disposal 
+ Printing labels 
+ Washing glassware


## Emergency contacts

Will Shoemaker 

+ Phone: 703-554-9264
+ E-mail: wrshoema@umail.iu.edu

Dr. Jay T. Lennon

+ lennonj@indiana.edu 