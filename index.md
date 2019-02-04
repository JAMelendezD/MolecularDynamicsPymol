# MNIST

This is a trained network to distinguish written digets from 0 to 9 from any image. The notebook loads the model (99.17% accuracy) trained with the MNIST dataset. Scripts with the model and the training can be found. Three example images are provided to test the network, the images or not from the trainig set nor the testing set. To test it images have to be written on a black background with a white brush.

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

