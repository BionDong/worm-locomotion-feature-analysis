/*
* Copyright 2019 Xianke Dong et al <xianke.dong@mail.mcgill.ca>
* This file is the image processing algorithm of the RoboWorm system.
* Specifically, it is developed to analyse the length/width/volume/
* curvature of a free moving C. elegans, and target an interested
* part on the worm body with high speed.
*
* This algorithm is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY for a perticular use. For the most up to
* date version of this software, see:
* https://github.com/BionDong/worm-locomotion-feature-analysis
*
* NOTE: If you use any portion of this code in your research, kindly cite:
* [1] Xianke Dong, Pengfei Song, and Xinyu Liu. "An Automated 
* Microfluidic System for Morphological Measurement and Size-Based 
* Sorting of C. Elegans." IEEE transactions on nanobioscience (2019),
*
* and our following works on C. elegans light-driven microrobot:
* [2] Xianke Dong, Sina Kheiri, Yangning Lu, Zhaoyi Xu, Mei Zhen, 
* Xinyu Liu, "Toward a living soft microrobot through optogenetic 
* locomotion control of Caenorhabditis elegans." Science Robotics 6, 
* eabe3950 (2021).

*/

/*
*  Created on: Oct. 10, 2019
*      Author: Xianke (Bion) Dong
*/

/*
* This project is built based on OPENCV with version of 3.0.0.
* The core algorithms for the morphorlogic feature analysis of free
* moving C. elegans are enclosed in the class named "WormTrack". In
* this code, the usesage of the WormTrack is illustrated by analysing
* 49 images for continuous worm moving.
*/
