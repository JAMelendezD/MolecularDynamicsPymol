# Prerequisites

The only assumption we make for this tutorial is that you have GROMACS and Pymol installed in your computer. For our case we are using Ubuntu 18.04, GROMACS 2018.3 and  Pymol 1.8.

# Pymol Basics

## Protein

We are gonna start with some pymol basics based on a enzyme called Mouse Thymidylate Synthase this enzyme catalyzes the conversion of deoxyuridine monophosphate (dUMP) to deoxythymidine monophosphate (dTMP). So the first thing we are gonna do is start a log file, download the protein from the protein data bank, remove the solvent and visualize it. Using the pymol terminal type: 

```
log_open log.pml
fetch 6F6Z
hide
remove solvent
show cartoon, 6F6Z
```
Most of the features in pymol can be typed in or selected in the Pymol viwer. Say for example we would like to change the color of our protein we could type: color cyan, 6F6Z or we could select in the Viwer by going to 6F6Z then clicking c (color) going down cyans and then cyan. For now we are gonna select 'c' then 'by ss' and finally 'Helix (cyan) Sheet (magenta) Loop (pink)'. Notice that if you select color by chain you can see that our protein is composed of two chains.

<p align="center">
  <img width="800" src="./media/color.png">
</p>

As you can see now we have our protein not only showing its secondary structure but also in different color so we can distinguish it. Now we are gonna focus on the ligands click on the letter S of the pymol viwer at the bottom right side of the screen a scroll bar should appear on top of the screen, her you can look and select every residue of your .pdb file. Scroll all the way to the right, you should see 4 particular residues named NOH and TGQ. 


<p align="center">
  <img width="800" src="./media/residues.png">
</p>

Now that we know the names of our ligands we are gonna select them using the terminal.

```
select lig1, resn NOH
select lig2, resn TGQ
```
You should know see in the viwer two new rows under 6F6Z, notice that know we can apply individal actions to our new selections. In this case we used resn or resname to select our residues we could have also used resi wich point towards the specific residue number. Now we are gonna show our ligands as spheres. 

```
show spheres, lig1
show spheres, lig2
set sphere_scale, 0.7
```

And color them using the viwer selecting color by element and select any color we like. Notice we also changed the size of our spheres from a default of 1 to 0.7. Finally we are gonna select our protein (6F6Z) an show it acessible surface. In the viwer select the letter s (show) on the top right and select surface.  

<p align="center">
  <img width="800" src="./media/surface.png">
</p>

Now we are ready to create a nice picture. Change the color of the protein to gray80 here, since we never defined the protein properly we have to change the color of our ligands again to any color by element. Now we are gonna look at the options we have on top of the pymol terminal window. Select display then background, set the background color to white and unselect the opaque option. Now select display then quality and select max quality. Now go to setting then transparency then surface and select 50. Finally go to scene then cache and select optimize. Now we are gonna ray trace our image, this setting has many options depending on what you want so we suggest reading the documentation in pymol. For know lest only test ray_trace_mode.

```
set ray_trace_mode, 1
ray 2000
```
Now that our image is ray traced we cannot move it otherwise we would have to repeat the last line of code to save the image just type 'png filename.png' now that you saved it you can move around to get another snap not forgetting to ray trace before saving. Notice that we used ray 2000 this creates an image of 2000x1356 in this case we only controlled the width but you can control both dimensions if you type 'ray 2000,1356' change the ray_trace_mode from 0 to 3 an see what type of images you can get.


<p align="center">
  <img width="800" src="./media/example2.png">
</p>

<p align="center">
  <img width="800" src="./media/example1.png">
</p>


# g_lomepro

The first thing we are going to do is download g_lomepro<sup>[1](#footnote1)</sup> a software developed by Vytautas Gapsys, Bert L. de Groot, Rodolfo Briones to calculate local properties of membranes. This link will take you to their website: [g_lomepro] (http://www3.mpibpc.mpg.de/groups/de_groot/g_lomepro.html). Once you downloaded it unzip it and  


<a name="footnote1">1</a>: Check the article about g_lomepro. [Vytautas Paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3882000/ "ncbi"). 

