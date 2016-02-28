% manual layer info:
% based on: /home/jason/projects/aro/ACd-data/for jason 151014 - with layer boundaries.pdf
L23cells=[103 69 89 88 108 112 58 116 86 91 59 104 65 56 66 97 111 76];
L3bcells=[176 75 77 80 172 171 167 170 114 83];
L5acells=[105 113 81 82 177 95];
L5bcells=[64 85 168 57 87 62 165];
L6cells=[62.5 60 67 61 166 72 117 74 71 73 55 102 109 169 84 68 70 96 94];

% Define 3 "layers" to process (super,mid,deep)
super=L23cells; 
mid=[L3bcells L5acells L5bcells]; 
deep=L6cells;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEW CLASSES BASED ON 5 IPs ACROSS 61 CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
  based on: 
    {'AHP(time2trough)','RMP','Ih','median(ISIs)','ISI2/ISI3'} from
    IP file: '/home/jason/models/dnsim/Tallie-ACd-cells/distributions/tallie_distributions_1.mat'
%}
% 5 clusters (100% cells clustered): IPs (_AHPtime2trough_RMP_medianISIs_Ih_ISI2-ISI3) (61 cells)
clusters{1} = [56 57 58 59 61 64 65 67 68 69 73 75 80 89 104 111 112 114 166 167 176 ];
clusters{2} = [55 60 70 74 85 88 92 94 97 108 109 168 177 ];
clusters{3} = [62 76 77 81 96 105 113 ];
clusters{4} = [66 82 83 84 86 117 ];
clusters{5} = [71 72 87 91 95 101 102 103 116 165 169 170 171 172 ];

% --------------------------------------
% MODEL CELLS:
%   Base on classes 1, 2, and 4
% MODEL LAYERS:
%   Superficial: Class 1 dominates, with a little class 4
%   Deep layer: Classes 1, 2, and 4 contribute equally
% --------------------------------------

% Class 1 has long AHP, slow adaptation, and dominates superficial layers.
% Class 2 has short AHP, rapid adaptation, and dominates deep layers.
% Class 3 has longest AHP, largest Ih sag, 
% Class 4 is rare and evenly distributed, with fastest/no AHP, elevated RMP
%         and abrupt stopping (few spikes, with ADP)
% Class 5 may have bigger Ih and slightly faster adaptation than Class 1, 
%         but otherwise similar IPs; is more concentrated in middle layers.

% Potential relation to Fiona's manual Groups:
% Group 1 <-> Class 3
% Group 2 <-> Class 1 (and maybe Class 5 in middle layers)
% Group 3 <-> Class 2
% Group 4 <-> Class 4

%{
Intrinsic Properties by Cell Class (mean/std,N)
                AHP(time2trough)	SpikeAmp        SpikeWidth        ThreshRate      RMP             1/ISI1          median(ISIs)		Ih                Ih-decay        ISI2/ISI3		
Class1 (n=21): 	(31  /12 ,21)			(66  /12 ,21)		(1.5 /0.28,21)		(1.8 /1.3,21)		(-75 /10 ,21)		(32  /13 ,21)		(11  /3.8,19)		(0.88/0.66,14)		(85  /63 ,14)   (0.9 /0.21,21)	
Class2 (n=13): 	(21  /12 ,13)			(68  /12 ,13)		(1.5 /0.32,13)		(1.7 /2  ,11)		(-75 /6.4,13)		(27  /15 ,13)		(22  /8.8,12)		(0.95/0.44,12)		(110 /67 ,12)		(0.56/0.13, 9)	
Class3 (n=7) : 	(51  /20 , 7)			(66  /8.4, 7)		(1.4 /0.17, 7)		(3.6 /1.5, 7)		(-71 /5.9, 7)		(37  /24 , 7)		(12  /7.3, 5)		(3.4 /1.1 , 7)		(130 /40 , 7)		(0.92/0.2 , 6)	
Class4 (n=6) : 	(9.8 /7.6, 3)			(59  /4.4, 3)		(1.5 /0.42, 3)		(2.6 /3.3, 3)		(-64 /3.3, 6)		(26  /17 , 6)		(31  /9.1, 6)		(2.3 /0.59, 6)		(92  /53 , 6)		(NaN /NaN , 0)	
Class5 (n=14): 	(39  /20 ,14)			(69  /4.4,14)		(1.4 /0.24,14)		(1.4 /1  ,13)		(NaN /NaN, 0)		(33  /29 ,14)		(14  /10 ,14)		(1.3 /0.95, 8)		(110 /48 , 8)		(0.73/0.18,11)	

################################################################
Laminar composition (within-layer percent of cells from each class):
          sup   mid   deep
Class1:   50%   30%   26%	
Class2:   17%   13%   32%	
Class3:   6%    22%   5%	
Class4:   11%   9%    11%	
Class5:   17%   26%   21%	
--------------------------
TotalNum: 18    23    19  (# cells per layer)
################################################################
Distribution across layers (for each cell class):
      Class1	Class2	Class3	Class4	Class5
sup :	43%     23%     14%     33%     21%	
mid :	33%     23%     71%     33%     43%	
deep:	24%     46%     14%     33%     29%	
--------------------------------------------
Tot#:	21      13      7       6       14    (# cells per class)
################################################################

Significant differences between cell classes:
#12-(Class1,Class2):
          ISI2/ISI3: p=2.4e-05
       median(ISIs): p=0.0019
   AHP(time2trough): p=0.034
#13-(Class1,Class3):
                 Ih: p=0.0003
         ThreshRate: p=0.019
         AHP(0-5ms): p=0.028
   AHP(time2trough): p=0.044
        AHP(5-20ms): p=0.046
#14-(Class1,Class4):
                RMP: p=0.00027
                 Ih: p=0.00079
       median(ISIs): p=0.003
   AHP(time2trough): p=0.018
#15-(Class1,Class5):
          ISI2/ISI3: p=0.025
#23-(Class2,Class3):
                 Ih: p=0.0005
          ISI2/ISI3: p=0.0047
   AHP(time2trough): p=0.0083
         AHP(0-5ms): p=0.018
         ThreshRate: p=0.034
       median(ISIs): p=0.045
#24-(Class2,Class4):
                RMP: p=9.5e-05
                 Ih: p=0.0013
#25-(Class2,Class5):
   AHP(time2trough): p=0.012
          ISI2/ISI3: p=0.026
#34-(Class3,Class4):
   AHP(time2trough): p=0.0018
       median(ISIs): p=0.0047
                RMP: p=0.018
                 Ih: p=0.032
#35-(Class3,Class5):
                 Ih: p=0.0017
         ThreshRate: p=0.0067
#45-(Class4,Class5):
   AHP(time2trough): p=0.0024
       median(ISIs): p=0.0052
           SpikeAmp: p=0.037
                 Ih: p=0.046


Intrinsic Properties (mean ? std) sorted by Heterogeneity (=std/mean)
Class1 (n=21):
(0.13):              RMP: -74.9 +/- 9.96            (n=21)
(0.18):         SpikeAmp: 65.6 +/- 11.7            (n=21)
(0.19):       SpikeWidth: 1.49 +/- 0.276           (n=21)
(0.22):      AHP(5-20ms): 1.01e+05 +/- 2.28e+04        (n=21)
(0.23):        ISI2/ISI3: 0.899 +/- 0.209           (n=21)
(0.33):     median(ISIs): 11.5 +/- 3.78            (n=19)
(0.39): AHP(time2trough): 30.7 +/- 12              (n=21)
(0.41):           1/ISI1: 31.9 +/- 12.9            (n=21)
(0.72):       ThreshRate: 1.82 +/- 1.31            (n=21)
(0.74):         Ih-decay:   85 +/- 62.9            (n=14)
(0.75):               Ih: 0.885 +/- 0.66            (n=14)
Class2 (n=13):
(0.085):              RMP: -75.4 +/- 6.44            (n=13)
(0.17):         SpikeAmp: 67.7 +/- 11.8            (n=13)
(0.22):       SpikeWidth: 1.49 +/- 0.321           (n=13)
(0.24):        ISI2/ISI3: 0.565 +/- 0.134           (n=9)
(0.28):      AHP(5-20ms): 1.02e+05 +/- 2.81e+04        (n=13)
(0.4):     median(ISIs): 21.7 +/- 8.76            (n=12)
(0.46):               Ih: 0.95 +/- 0.437           (n=12)
(0.54): AHP(time2trough): 21.4 +/- 11.7            (n=13)
(0.58):           1/ISI1: 26.7 +/- 15.4            (n=13)
(0.59):         Ih-decay:  113 +/- 66.5            (n=12)
(1.2):       ThreshRate: 1.69 +/- 1.99            (n=11)
Class3 (n=7):
(0.082):              RMP: -71.3 +/- 5.87            (n=7)
(0.12):       SpikeWidth: 1.41 +/- 0.166           (n=7)
(0.13):         SpikeAmp: 66.5 +/- 8.36            (n=7)
(0.21):        ISI2/ISI3: 0.923 +/- 0.198           (n=6)
(0.31):               Ih: 3.43 +/- 1.05            (n=7)
(0.31):         Ih-decay:  129 +/- 40.2            (n=7)
(0.33):      AHP(5-20ms): 7.66e+04 +/- 2.55e+04        (n=7)
(0.41): AHP(time2trough): 50.5 +/- 20.5            (n=7)
(0.41):       ThreshRate: 3.61 +/- 1.49            (n=7)
(0.6):     median(ISIs): 12.1 +/- 7.29            (n=5)
(0.66):           1/ISI1: 36.6 +/- 24.4            (n=7)
Class4 (n=6):
(0.052):              RMP: -63.9 +/- 3.34            (n=6)
(0.074):         SpikeAmp: 58.8 +/- 4.37            (n=3)
(0.26):               Ih: 2.27 +/- 0.59            (n=6)
(0.28):       SpikeWidth: 1.53 +/- 0.421           (n=3)
(0.3):     median(ISIs): 30.5 +/- 9.08            (n=6)
(0.5):      AHP(5-20ms): 9.37e+04 +/- 4.71e+04        (n=3)
(0.58):         Ih-decay: 92.2 +/- 53.5            (n=6)
(0.67):           1/ISI1: 25.7 +/- 17.3            (n=6)
(0.78): AHP(time2trough):  9.8 +/- 7.64            (n=3)
(1.3):       ThreshRate: 2.57 +/- 3.28            (n=3)
(NaN):        ISI2/ISI3:  NaN +/- NaN             (n=0)
Class5 (n=14):
(0.064):         SpikeAmp: 68.8 +/- 4.43            (n=14)
(0.18):       SpikeWidth: 1.39 +/- 0.244           (n=14)
(0.24):        ISI2/ISI3: 0.733 +/- 0.176           (n=11)
(0.28):      AHP(5-20ms): 8.64e+04 +/- 2.43e+04        (n=14)
(0.44):         Ih-decay:  109 +/- 48.4            (n=8)
(0.51): AHP(time2trough): 38.6 +/- 19.8            (n=14)
(0.69):     median(ISIs): 14.5 +/- 9.98            (n=14)
(0.7) :               Ih: 1.35 +/- 0.948           (n=8)
(0.71):       ThreshRate: 1.41 +/- 1               (n=13)
(0.89):           1/ISI1: 32.7 +/- 29.2            (n=14)
(NaN):              RMP:  NaN +/- NaN             (n=0)

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ORIGINAL CLASSES BASED ON 7 IPs ACROSS 61 CELLS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
  based on: {'AHP(total)','AHP(time2trough)','SpikeWidth','ThreshRate','RMP','1/ISI1','ISI2/ISI3'}
    i.e.,
      'AHP total'
      'AHP time to trough'
      'Spike width at half height'
      'Spike rate at threshold'
      'RMP'
      'Initial spike instant frequ...'
      'Ratio of 2nd to 3rd ISI'
    IP file: '/home/jason/projects/aro/analysis/clustering/tallie_IPs_150410.mat'
%            '/home/jason/projects/aro/analysis/clustering/tallie_IP_labels_150420.mat'
%}  
% cell IDs in each class
clusters{1} = [56 59 60 61 73 92 97 101 103 109 111 112 114 165 166 168 170 171 172 177 ];
clusters{2} = [55 57 58 62 65 67 68 69 70 71 74 80 83 84 86 88 94 96 104 167 169 ];
clusters{3} = [72 75 87 95 102 105 113 116 176 ];
clusters{4} = [76 77 81 89 108 117 ];
clusters{5} = [64 66 82 85 91 ];
%{
Intrinsic Properties by Cell Class (mean/std,N)
          AHP(time2trough)	SpikeAmp        SpikeWidth      ThreshRate        RMP             1/ISI1          median(ISIs)		Ih                Ih-decay          ISI2/ISI3		
Class1: 	(28  /8.9,20)			(72  /10 ,20)		(1.4 /0.2 ,20)	(1.3 /0.79,18)		(-81 /5.3,14)		(36  /25 ,20)		(18  /11 ,18)		(0.68/0.28,12)		(96  /71 ,12)     (1.6 /3.7  ,18)	
Class2: 	(23  /11 ,19)			(64  /7.7,19)		(1.4 /0.18,19)	(1.6 /0.82,19)		(-69 /4.8,19)		(28  /17 ,21)		(16  /8.9,21)		(1.5 /1.1 ,19)    (101 /54 ,19)     (0.83/0.26 ,18)	
Class3: 	(62  /11 , 9)			(69  /6.9, 9)		(1.5 /0.24, 9)	(1.8 /1.1 , 9)    (-80 /5.3, 4)		(26  /17 , 9)		(8.8 /2.9, 8)		(1.6 /0.81, 7)		(126 /53 , 7)     (0.8 /0.085, 7)	
Class4: 	(34  /23 , 6)			(59  /9.9, 6)		(1.5 /0.3 , 6)	(5.9 /0.98, 6)		(-70 /9.8, 6)		(40  /17 , 6)		(18  /10 , 4)		(4.1 /0.63, 4)		(87  /36 , 4)     (0.96/0.22 , 4)	
Class5: 	(24  /17 , 4)			(59  /6.3, 4)		(2.1 /0.19, 4)	(0.84/0.98, 3)		(-64 /6.2, 4)		(21  /9.3, 5)		(26  /8.7, 5)		(1.6 /0.44, 5)		(113 /65 , 5)     (0.77/0    , 1)	

          AHP(0-5ms)              AHP(5-20ms)             AHP(20-200ms)         AHP(total)              AHP(trough)		
Class1: 	(2.2e+04/5.2e+03,20)		(9.9e+04/2.2e+04,20)		(9.6e+05/3e+05,20)		(1.1e+06/3.3e+05,20)		(8.4 /2.3,20)	
Class2: 	(2.7e+04/5.9e+03,19)		(1.1e+05/2.1e+04,19)		(9.4e+05/2.1e+05,19)  (1.1e+06/2.3e+05,19)		(8.8 /2.4,19)	
Class3: 	(1.5e+04/5.7e+03, 9)		(7.3e+04/2.1e+04, 9)		(1.1e+06/1.4e+05, 9)	(1.2e+06/1.6e+05, 9)		(8.9 /1.7, 9)	
Class4: 	(1.7e+04/9.6e+03, 6)		(6.7e+04/2.9e+04, 6)		(7e+05/3.1e+05, 6)		(7.9e+05/3.5e+05, 6)		(4.9 /2.7, 6)	
Class5: 	(2.3e+04/4.7e+03, 4)		(9.4e+04/3.2e+04, 4)		(9.3e+05/4e+05, 4)		(1e+06/4.4e+05, 4)    	(7.4 /2.2, 4)	

################################################################
Laminar composition (within-layer percent of cells from each class):
          sup   mid   deep
Class1:   33%   30%   26%	
Class2:   33%   22%   53%	
Class3:   6%    26%   11%	
Class4:   17%   9%    5%	
Class5:   11%   13%   0%	
--------------------------
TotalNum: 18    23    19  (# cells per layer)
################################################################
Distribution across layers (for each cell class):
        Class1	Class2	Class3	Class4	Class5
sup :   30%     29%     11%     50%     40%	
mid :   35%     24%     67%     33%     60%	
deep:   25%     48%     22%     17%     0%	
--------------------------------------------
Tot#:   20      21      9       6       5   (# cells per class)
################################################################

Significant differences between cell classes:
(Class1,Class2):
              RMP: p=1.3e-07
               Ih: p=0.0042
         SpikeAmp: p=0.0061
       AHP(0-5ms): p=0.0079
(Class1,Class3):
  AHP(time2trough): p=1.7e-06
     median(ISIs): p=0.0022
       AHP(0-5ms): p=0.0081
      AHP(5-20ms): p=0.0089
               Ih: p=0.021
(Class1,Class4):
       ThreshRate: p=1.1e-05
               Ih: p=0.0011
      AHP(trough): p=0.024
         SpikeAmp: p=0.022
              RMP: p=0.032
      AHP(5-20ms): p=0.046
(Class1,Class5):
       SpikeWidth: p=0.0013
               Ih: p=0.0052
              RMP: p=0.0056
         SpikeAmp: p=0.011
           1/ISI1: p=0.036
(Class2,Class3):
  AHP(time2trough): p=1.4e-07
       AHP(0-5ms): p=0.00011
      AHP(5-20ms): p=0.00069
     median(ISIs): p=0.0045
    AHP(20-200ms): p=0.0088
              RMP: p=0.016
(Class2,Class4):
       ThreshRate: p=2e-05
               Ih: p=0.00032
      AHP(5-20ms): p=0.015
      AHP(trough): p=0.016
(Class2,Class5):
       SpikeWidth: p=0.002
(Class3,Class4):
       ThreshRate: p=6.8e-06
               Ih: p=0.00059
      AHP(trough): p=0.014
    AHP(20-200ms): p=0.019
  AHP(time2trough): p=0.029
(Class3,Class5):
       SpikeWidth: p=0.0017
              RMP: p=0.0082
     median(ISIs): p=0.01
  AHP(time2trough): p=0.015
       AHP(0-5ms): p=0.032
         SpikeAmp: p=0.038
(Class4,Class5):
               Ih: p=0.0011
       ThreshRate: p=0.0017
       SpikeWidth: p=0.0034
           1/ISI1: p=0.039


Intrinsic Properties (mean ? std) sorted by Heterogeneity (=std/mean)
Class1 (n=20):
(0.0653):             RMP: -81.4 ? 5.32           (n=14)
(0.141):         SpikeAmp: 72.1 ? 10.2            (n=20)
(0.146):       SpikeWidth: 1.35 ? 0.198           (n=20)
(0.222):      AHP(5-20ms): 9.86e+04 ? 2.19e+04    (n=20)
(0.315): AHP(time2trough): 28.3 ? 8.91            (n=20)
(0.42):                Ih: 0.676 ? 0.284          (n=12)
(0.587):     median(ISIs): 18.2 ? 10.7            (n=18)
(0.618):       ThreshRate: 1.29 ? 0.794           (n=18)
(0.681):           1/ISI1: 36.2 ? 24.6            (n=20)
(0.737):         Ih-decay: 96.3 ? 71              (n=12)
(2.32):         ISI2/ISI3: 1.61 ? 3.74            (n=18)
Class2 (n=21):
(0.0702):             RMP: -68.6 ? 4.81           (n=19)
(0.121):         SpikeAmp: 63.7 ? 7.7             (n=19)
(0.125):       SpikeWidth: 1.41 ? 0.177           (n=19)
(0.194):      AHP(5-20ms): 1.09e+05 ? 2.11e+04    (n=19)
(0.32):         ISI2/ISI3: 0.826 ? 0.265          (n=18)
(0.472): AHP(time2trough): 22.6 ? 10.7            (n=19)
(0.497):       ThreshRate: 1.64 ? 0.816           (n=19)
(0.537):         Ih-decay:  101 ? 54.4            (n=19)
(0.57):      median(ISIs): 15.6 ? 8.86            (n=21)
(0.592):           1/ISI1: 27.9 ? 16.5            (n=21)
(0.707):               Ih: 1.49 ? 1.05            (n=19)
Class3 (n=9):
(0.0663):             RMP:  -80 ? 5.31            (n=4)
(0.1):           SpikeAmp: 68.8 ? 6.91            (n=9)
(0.106):        ISI2/ISI3: 0.803 ? 0.0848         (n=7)
(0.156):       SpikeWidth: 1.51 ? 0.235           (n=9)
(0.174): AHP(time2trough): 61.7 ? 10.8            (n=9)
(0.282):      AHP(5-20ms): 7.35e+04 ? 2.07e+04    (n=9)
(0.327):     median(ISIs): 8.78 ? 2.87            (n=8)
(0.417):         Ih-decay:  127 ? 52.8            (n=7)
(0.499):               Ih: 1.62 ? 0.809           (n=7)
(0.604):       ThreshRate: 1.77 ? 1.07            (n=9)
(0.652):           1/ISI1: 25.8 ? 16.8            (n=9)
Class4 (n=6):
(0.141):              RMP: -69.7 ? 9.82           (n=6)
(0.156):               Ih: 4.05 ? 0.631           (n=4)
(0.166):       ThreshRate:  5.9 ? 0.978           (n=6)
(0.167):         SpikeAmp: 59.2 ? 9.88            (n=6)
(0.205):       SpikeWidth: 1.46 ? 0.299           (n=6)
(0.227):        ISI2/ISI3: 0.963 ? 0.219          (n=4)
(0.412):           1/ISI1:   40 ? 16.5            (n=6)
(0.416):         Ih-decay: 86.8 ? 36.1            (n=4)
(0.438):      AHP(5-20ms): 6.7e+04 ? 2.93e+04     (n=6)
(0.564):     median(ISIs):   18 ? 10.2            (n=4)
(0.683): AHP(time2trough): 33.6 ? 22.9            (n=6)
Class5 (n=5):
(  0):          ISI2/ISI3: 0.772 ? 0              (n=1)
(0.0886):      SpikeWidth: 2.09 ? 0.185           (n=4)
(0.0969):             RMP:  -64 ? 6.2             (n=4)
(0.108):         SpikeAmp: 58.6 ? 6.33            (n=4)
(0.27):                Ih: 1.64 ? 0.444           (n=5)
(0.339):     median(ISIs): 25.6 ? 8.7             (n=5)
(0.341):      AHP(5-20ms): 9.4e+04 ? 3.21e+04     (n=4)
(0.454):           1/ISI1: 20.6 ? 9.34            (n=5)
(0.571):         Ih-decay:  113 ? 64.5            (n=5)
(0.726): AHP(time2trough): 23.8 ? 17.2            (n=4)
(1.16):        ThreshRate: 0.841 ? 0.979          (n=3)

%}