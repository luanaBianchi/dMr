# dMr
- Enter the name of an input data file (two columns, M vs H). It can contain 
a single recoil loop traced 
after saturation (either positive or negative), or 
a sequence of such loops. In case 
a major loop is also included, it should 
be the last one in the data sequence.


- The program offers two options:
1) The data file contains recoil loop(s) only. In case, the shift values along 
the H- and M- axes (if any) must be provided.
2) The data file represents a sequence of recoil loop(s) together with the 
major hysteresis loop. In this case, an image will pop-up showing the major 
loop, the coercive and swithching fileds and the corresponding shift values. 
The code will suggest the M- and H-field (Heb) shifts obtained from the 
coercive fields; however, other shift values can be used.

- The next line is self-sufficient.

- The dMr's field step can be changed, though the value proposed 
usually is good enough.


- A graph will pop up in a new window. It can be saved as an image. 

- 

Four data files will be generated in the same directory:
 
1) centered_loop(s).dat  - the centralized data: H x M; 

2) dMrs.dat - it contains a single dMr plot or the sequence of such 
    plots if the 
input data file contains a sequence of loops. The 3rd 
    column gives the recoil-field value corresponding to each dMr plot; 

3) major-loop_parameters.dat  - the values of the swithching and 
    coercive fields, the M- and H-shifts, the remannt magnetizations, 
    the recoil field(s) and the extreme values (with the respective fields) 
    of each dMr plot. 
4) area.dat - there are seven columns in this file, i.e.,  
    - the value of the recoil field;
    - the value of the field corresponing to the maximum of dMr(H);
    - the value of the magnetization, Mmax, at the above field;
    - the value of the field corresponing to the minimum of dMr(H);
    
- the value of the magnetization, Mmin, at the above field;

    - the value of the area with positive M-values;

    - the value of the area with negative M-values;

In the compressed file there are two sample images corresponding 
to the results the code gives 
when the two input sample files 
("recoil-loop.dat" and "recoil-&-major-loops.dat") are used. The first 
file contains one recoil loop only, and the second one contains a 
recoil loop together with 
the respective major loop (in this example, 
the major loop is asymmetric and shifted by -215 Oe).



We welcome suggestions for bug reports and improvements. 

Please, contact Julian Geshev (julian@if.ufrgs.br). 
