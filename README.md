There are two parts which I implement parallelization. The first one is computeForce() function.
Every thread owns its indexStart and indexEnd.  When multiple threads run in the function, the data will only be modified by one thread.
The second one is applyForce().
Same as computingForce, the data won’t be modified by multiple threads because that each thread own its working section. 


# Partition Idea
The partition idea comes from the HW1. Basically, I fairly allocate threads to get same amount of data. However, if there remain some, I will give every thread another one.
This is the function that I divide data into threads.

# Synchronization Problem
Each thread has the following works: buildingTree(), updateNodePos(), computeForce(), applyForce(), updateBoundWidth() and deleteTree(). The parallel functions are applyForce() and computeForce(), and others are sequential functions. I will add barrier after all function for secure all threads finish and go to next function. But when it is the sequential one, I will add mutex and flag variable to limit only one thread can enter the function.
For example, the above one is the sequential buildingTree() function. I use the flag variable buildingTreeFlag and check whether it is i or not. If it is, add one, that can prevent other thread in the i step enter the same function again.


# Technique for Reducing Execution Time and increasing Scalability
I will handle the real position for each node in the tree at last. The updateNodePos() will let me traverse the node and handle the position correctly.
Tree should be built every step. Root width will change step by step. The applyForce() not only help me to update the position but also collect every threads’ bound width. Though the updateBoundWidth() is sequential function, but traversing the amount of thread is much faster than traversing whole leaf in the node.

