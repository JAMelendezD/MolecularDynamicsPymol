# Prerequisites

The only assumption we make for this tutorial is that you have GROMACS and Pymol installed in your computer. For our case we are using Ubuntu 18.04, GROMACS 2018.3 and  Pymol 1.8.

# Pymol Basics

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
Now that our image is ray traced we cannot move it otherwise we would have to repeat the last line of code to save the image just type 'png filename.png' now that you saved it you can move around to get another snap not forgetting to ray trace before saving. Notice that we used ray 2000 this creates an image of 2000x 


# g_lomepro

The first thing we are going to do is download g_lomepro<sup>[1](#footnote1)</sup> a software developed by Vytautas Gapsys, Bert L. de Groot, Rodolfo Briones to calculate local properties of membranes. This link will take you to their website: [g_lomepro] (http://www3.mpibpc.mpg.de/groups/de_groot/g_lomepro.html). Once you downloaded it unzip it and  



```
if (isAwesome){
  return true
}
```

<p align="center">
  <img width="800" src="./media/mnist.png">
</p>

# Dogs and Cats

This is a trained network to distinguish between cats and dogs from any image. The notebook loads the model (85% accuracy) trained with the kaggle dataset. Scripts with the actual network and data pre-processing can also be found. Two example images are provided to test the network, the images or not from the trainig set nor the validation set.

<p align="center">
  <img width="200" height="200" src="./media/dog.jpg">
  <img width="200" height="200" src="./media/cat.jpg">
</p>


# Visualization

So far we can see that for training on an image database we normally use a deep convolutional neural network were we can also conbine dense layers at the end before the prediction. This method has provided a very effective way to train neural networks, in fact there exists pre-trained nets based on this idea for example the VGG16<sup>[1](#footnote1)</sup>.

<p align="center">
  <img width="800" src="./media/vgg16.png">
</p>

The idea behind convolutional neural networks sometimes is hard to picture so the purpose of this scripts is to get an introduction into images and openCV to then visualize the process behind convolutional neural networks. We start with our previously trained model of cats and dogs but this time our objective is the intermediate steps not the final prediction. We know convolutional layers are based on filters.

<p align="center">
  <img width="400" src="./media/filters.png">
</p>

Filters are just an array of numbers that we can visualzie with colors in the same way we do with an image, this set of numbers are optimize during the training of the network. We apply the filters to the images to extract particular characteristics for example applying the previous filter to an image of a dog we obtain:

<p align="center">
  <img width="400" src="./media/featuremap.png">
</p>

Now instead of one image we have a number of images equal to the number of filters these new images are normally refered as feature-maps. We can also visualize this as being a 3-dimensional set of data one for the dog (left) and one for the cat (right) we can already begin to see the diferences of a cat an a dog to the eyes of the network.

<p align="center">
  <img width="800" src="./media/3dout.png">
</p>

We can then keep applying convolutional layers to focus in different characteristics that are unique to each feature map. This process is applied multiple times with intermidiate pooling layers that decrease the size of the outputs until eventually we are left with a particular image size that we want to flatten to then pass it through a dense layer. This dense layer consist on weights and biases in our particular case this were our dense layer weights:

<p align="center">
  <img width="600" src="./media/weights.png">
</p>

Now the network is ready to give an answer it takes those final weights and multiplies them with the input. Here is the diferences of a cat and a dog in this final step before the weights and after the weights. The blue dots represent the dog and the red ones the cat. 

<p align="center">
  <img width="800" src="./media/differences.png">
</p>

Finally since it is a binary network we apply the sigmoid function to get our particular result. If we get a result greater than 0.5 the prediction is a cat if it is less that 0.5 then we get a prediction of a dog. 

<p align="center">
  <img width="700" src="./media/results.png">
</p>

We can see the importance of the weights and the sigmoid function in our neural network, without this, the distinction of the animals from just the convolutional output would be impossible.

# Raspberry Pi  

When working with neural networks there is always a question in the back of everyones mind, How deep can we go? The answer to that question depends on many aspects, in particular computation power and expected results.

<p align="center">
  <img width="800" src="./media/xkcd.png">
</p>

As we can see our cat-dog network is optimized for performing better results that the raspberry network this is simply because of computation power. Since for the raspberry proyect we expect to determine Persons, cats and dogs from live video using a logitech camera, due to limitations in the processing power we design a small network compared to the cat-dog to achieve higher framerates in the rapberry. The main purpose is to run camara.ipynb as a .py file in the raspberry Pi with a webcam, to get live video that classifies the given input into the categories. (It was trained with persons doing actions such as applauding, writing, etc).  

<a name="footnote1">1</a>:To learn how to do this graphs you can visit:[HowToPlot](https://github.com/JAMelendezD/HowToGraph "Julian's Repository") or check the article [Rougier Plos One Paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003833 "Plos One"). 

