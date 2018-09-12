CloudAcademyAlbedo
------------------

Getting Started 
===============

1) Download the repository 
````````````````````````````

.. code-block:: bash
	
	git clone https://github.com/natashabatalha/CloudAcademyAlbedo.git
	cd CloudAcademyAlbedo

2) Install some Python stuff, if you don't have these 
````````````````````````````````````````````````````

.. code-block:: bash 

	pip install pandas numpy bokeh astropy scipy cython

OR if you are conda user 

.. code-block:: bash

	conda install pandas numpy bokeh astropy scipy cython

3) Download zip files and put them in the right spot 
````````````````````````````````````````````````````
Download the two files `at this link <https://drive.google.com/drive/folders/1Helb2qJ1s_lJUIAJbZKhuUSnY8BHT_hD?usp=sharing>`_ by right clicking on each file to download the whole zip.

.. code-block:: bash

	mv freedman CloudAcademyAlbedo/AlbedoCode/
	mv opacities CloudAcademyAlbedo/AlbedoCode/

4) COMPILE!
````````````

.. code-block:: bash 

	cd CloudAcademyAlbedo/AlbedoCode 
	make

5) Open up **FunWithAlbedoCode.ipynb** and follow the directions! 
`````````````````````````````````````````````````````````````````

.. code-block:: bash 

	cd CloudAcademyAlbedo
	jupyter notebook FunWithAlbedoCode.ipynb

