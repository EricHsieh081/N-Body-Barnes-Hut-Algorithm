#include <stdio.h>
#include <stdlib.h>
#include <X11/Xlib.h>
#include <unistd.h>
#include <cmath>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <assert.h>
#include <time.h>
using namespace std;

#define G 6.67384 
#define E 1.0E-6

int numOfStar;
int steps;
int numOfThread;
double mass;
double timeOfStep;
double GMMS;
double theta;
double xmin;
double ymin;
double coordiLength;
double xwindowLength;
double cor;
char *xwindowOpen;

struct star_Info {
	double x;
	double y;
	double vx;
	double vy;
};

struct star_Info *stars; // real data here
struct star_Info *tempstars; // temp place to place the new data

int *workingSection; // memorize the thread working section start place
int *workingStartPoint;

double boundMaxY, boundMinY, boundMaxX, boundMinX, boundWidth; // correct bound saving
double *boundMinXArr, *boundMaxXArr, *boundMinYArr, *boundMaxYArr; // it will be temp array

struct timespec start;

double gettime() {
	struct timespec now;
	clock_gettime(CLOCK_MONOTONIC, &now);
	double sec = now.tv_sec - start.tv_sec;
	double ns = now.tv_nsec - start.tv_nsec;

	return (sec * 1E3) + (ns * 1E-6);
}

struct bhTree {
	struct bhTree *parent;
	struct bhTree *child_NW;
	struct bhTree *child_NE;
	struct bhTree *child_SW;
	struct bhTree *child_SE;
	double nodeWidth;
	double massPosX;
	double massPosY;
	double centerX;
	double centerY;
	int starIndex;
	int starNumLower;
};

struct bhTree *root = NULL;
struct bhTree *cur = NULL;
struct bhTree *curParent = NULL;

GC gc;
Display *display;
Window window;      //initialization for a window
int screen;         //which screen 
const double PI = 3.14159;
double IOStart, buildTreeStart;
double IOTime, buildTreeTime, computingTime;

int updateNodePosFlag, updateBoundWidthFlag, deleteTreeFlag, buildTreeFlag, testFlag;
pthread_mutex_t mutex;
pthread_barrier_t barrier;

void initGraph(int width, int height);
void draw(int x, int y);
void allocateStar();
void *countingStar(void *threadID);
void updateNodePos();
void updateBoundWidth();
void deleteTree();
void applyForce(int indexStart, int indexEnd, long id);
void computeForce(int indexStart, int indexEnd);
void traverse(int starIndex, double &Fx, double &Fy, int forceORpos, bhTree *now, int level);
void buildTree(int indexStart, int indexEnd);
void create_node(int fall, int starIndex, bhTree **now, bhTree **nowParent);
int whereToFall(int starIndex, bhTree **now, bhTree **nowParent);

Pixmap buffer;

int main(int argc,char *argv[]) {

	IOTime = buildTreeTime = computingTime = 0;
	IOStart = buildTreeStart;
	clock_gettime(CLOCK_MONOTONIC, &start);

	numOfThread = atoi(argv[1]);
	mass = atof(argv[2]);
	steps = atoi(argv[3]);
	timeOfStep = atof(argv[4]);
	char *file = argv[5];
	theta = atof(argv[6]);
	xwindowOpen = argv[7];
	xmin = atof(argv[8]);
	ymin = atof(argv[9]);
	coordiLength = atof(argv[10]);
	xwindowLength = atof(argv[11]);

	cor = xwindowLength/coordiLength;

	GMMS = G*mass*timeOfStep*(1E-11); //G*M*M*timeOfStep/M simple

	/*read star and initialize array*/
	FILE *fp;
	fp = fopen(file, "r");
	if(fp == NULL) {
		perror("File doesn't exist.");
		fclose(fp);
		return 0;
	}
	fscanf(fp, "%d", &numOfStar);
	stars = (star_Info*)malloc(sizeof(star_Info)*numOfStar);
	tempstars = (star_Info*)malloc(sizeof(star_Info)*numOfStar);

	IOStart = gettime();
	fscanf(fp, "%lf %lf %lf %lf", &stars[0].x, &stars[0].y, &stars[0].vx, &stars[0].vy);
	IOTime = IOTime + gettime() - IOStart;
	boundMinX = stars[0].x;
	boundMinY = stars[0].y;
	boundMaxX = stars[0].x;
	boundMaxY = stars[0].y; 
	//printf("%lf %lf %lf %lf\n", stars[0].x, stars[0].y, stars[0].vx, stars[0].vy);
	boundWidth;
	for(int i = 1; i < numOfStar; i++) {
		IOStart = gettime();
		fscanf(fp, "%lf %lf %lf %lf", &stars[i].x, &stars[i].y, &stars[i].vx, &stars[i].vy);
		IOTime = IOTime + gettime() - IOStart;
		if(stars[i].x < boundMinX)
			boundMinX = stars[i].x;
		if(boundMaxX < stars[i].x)
			boundMaxX = stars[i].x;
		if(stars[i].y < boundMinY)
			boundMinY = stars[i].y;
		if(boundMaxY < stars[i].y)
			boundMaxY = stars[i].y;
		//printf("%lf %lf %lf %lf\n", stars[i].x, stars[i].y, stars[i].vx, stars[i].vy);
	}
	boundWidth = ((boundMaxX - boundMinX) > (boundMaxY - boundMinY)) ? (boundMaxX - boundMinX) : (boundMaxY - boundMinY);
	fclose(fp);

	/*handle how to divide star and balance thread*/
	allocateStar();

	/*initialize root & boundArr*/
	boundMaxYArr = (double*)malloc(sizeof(double)*numOfThread);
	boundMaxXArr = (double*)malloc(sizeof(double)*numOfThread);
	boundMinYArr = (double*)malloc(sizeof(double)*numOfThread);
	boundMinXArr = (double*)malloc(sizeof(double)*numOfThread);

	root = (bhTree *)malloc(sizeof(bhTree));
	root->parent = NULL;
	root->child_SE = NULL;
	root->child_SW = NULL;
	root->child_NE = NULL;
	root->child_NW = NULL;
	root->nodeWidth = boundWidth;
	root->starIndex = -1;
	root->starNumLower = 0;
	root->massPosX = 0;
	root->massPosY = 0;
	root->centerX = boundMinX + boundWidth/2;
	root->centerY = boundMinY + boundWidth/2;

	/*initalize pthread and wait until join*/
	if(!strcmp(xwindowOpen, "enable")) 	
		initGraph(xwindowLength, xwindowLength);

	if(numOfStar < numOfThread) 
		numOfThread = numOfStar;

	pthread_mutex_init(&mutex, NULL);
	pthread_barrier_init(&barrier, NULL, (unsigned) numOfThread);

	pthread_t threads[numOfThread];
	for(long id = 0; id < numOfThread; id++)
		pthread_create(&threads[id], NULL, countingStar, (void *)id);

    for (int i = 0; i < numOfThread; i++)
		pthread_join(threads[i], NULL);

	free(stars);
	free(tempstars);
	free(workingSection);
	free(workingStartPoint);

	pthread_mutex_destroy(&mutex);
	pthread_barrier_destroy(&barrier);
	computingTime = gettime() - buildTreeTime - IOTime;
	//printf("%g\n", gettime());
	printf("%g\t%g\t%g\n", computingTime, buildTreeTime, IOTime);
	return 0;
}

void allocateStar() {
	workingSection = (int*)malloc(sizeof(int)*numOfThread);
	workingStartPoint = (int*)malloc(sizeof(int)*numOfThread);
	int base = numOfStar/numOfThread;
	int left = numOfStar%numOfThread;

	for(int i = 0; i < numOfThread; i++) 
		workingSection[i] = base;
	for(int i = 0; i < left; i++) 
		workingSection[i]+=1;
	for(int i = 0; i < numOfThread; i++) 
		workingStartPoint[i] = 0;
	for(int i = 1; i < numOfThread; i++)
		workingStartPoint[i] = workingSection[i-1]+workingStartPoint[i-1];
}

void *countingStar(void *threadID) {
	long id = (long)threadID;
	int indexStart = workingStartPoint[id];
	int indexEnd = workingSection[id] + indexStart;

	//printf("I am Thread %2d, working on [%3d,%3d)\n", id, indexStart, indexEnd);

	for(int i = 0; i < steps; i++) {

		pthread_mutex_lock(&mutex);
		if(buildTreeFlag == i) {
			buildTreeFlag = i+1;
			buildTree(0, numOfStar); // parallel & for building tree
		}
		pthread_mutex_unlock(&mutex);
		pthread_barrier_wait(&barrier);

		pthread_mutex_lock(&mutex);
		if(updateNodePosFlag == i) {
			updateNodePosFlag = i+1;
			updateNodePos(); // sequential & for updating node position and draw line
		}
		pthread_mutex_unlock(&mutex);
		pthread_barrier_wait(&barrier);


		computeForce(indexStart, indexEnd); // parallel & for computing force to the tempstars
		pthread_barrier_wait(&barrier);

		applyForce(indexStart, indexEnd, id); //parallel & for apply Force back to the stars
		pthread_barrier_wait(&barrier);
		
		pthread_mutex_lock(&mutex);
		if(updateBoundWidthFlag == i) {
			updateBoundWidthFlag = i+1;
			updateBoundWidth(); //sequential & for updating boundWidth and draw star
		}
		pthread_mutex_unlock(&mutex);
		pthread_barrier_wait(&barrier); 
		
		pthread_mutex_lock(&mutex);
		if(deleteTreeFlag == i) {
			//printf("%ld - step: %d\n", id, i);
			deleteTreeFlag = i+1;
			deleteTree(); // sequential & for deleting tree, remain root
		}
		pthread_mutex_unlock(&mutex);
		pthread_barrier_wait(&barrier);
	}
	pthread_exit(NULL);
}

void updateNodePos() {
	double useless = 0;
	if(!strcmp(xwindowOpen, "enable")) {
		XSetForeground(display,gc,BlackPixel(display,screen));
		XFillRectangle(display,buffer,gc,0,0,xwindowLength,xwindowLength);
		XSetForeground(display,gc,0x606040);
	}
	traverse(-1, useless, useless, 1, root, 0);
}

void updateBoundWidth() {
	boundMaxX = boundMaxXArr[0];
	boundMaxY = boundMaxYArr[0];
	boundMinX = boundMinXArr[0];
	boundMinY = boundMinYArr[0];
	for(int i = 1; i < numOfThread; i++) {
		if(boundMaxXArr[i] > boundMaxX)
			boundMaxX = boundMaxXArr[i];
		if(boundMaxYArr[i] > boundMaxY)
			boundMaxY = boundMaxYArr[i];
		if(boundMinXArr[i] < boundMinX)
			boundMinX = boundMinXArr[i];
		if(boundMinYArr[i] < boundMinY)
			boundMinY = boundMinYArr[i];
	}
	boundWidth = ((boundMaxX - boundMinX) > (boundMaxY - boundMinY)) ? (boundMaxX - boundMinX) : (boundMaxY - boundMinY);
	if(!strcmp(xwindowOpen, "enable")) {
		IOStart = gettime();
		for(int j = 0; j < numOfStar; j++) {
			draw((int)((stars[j].x-xmin)*cor), (int)((stars[j].y-ymin)*cor));
		}
		XCopyArea(display, buffer, window, gc, 0, 0, xwindowLength, xwindowLength, 0, 0);
		XFlush(display);
		IOTime = IOTime + gettime() - IOStart;
	}			
}

void deleteTree() {
	struct bhTree *temp;
	int where;

	curParent = NULL;
	cur = root;
	while(1) {
		if(cur == root && cur->child_SE == NULL && cur->child_SW == NULL && cur->child_NE == NULL && cur->child_NW == NULL)
			break;
		if(cur->child_SE != NULL) {
			curParent = cur;
			cur = cur->child_SE;
			where = 4;
		}
		else if(cur->child_SW != NULL) {
			curParent = cur;
			cur = cur->child_SW;
			where = 3;
		}
		else if(cur->child_NW != NULL) {
			curParent = cur;
			cur = cur->child_NW;
			where = 2;
		}
		else if(cur->child_NE != NULL) {
			curParent = cur;
			cur = cur->child_NE;
			where = 1;
		}

		if(cur != root) {
			temp = cur;
			cur = curParent;
			curParent = cur->parent;
			temp->parent = NULL;
			free(temp);

			if(where == 1)
				cur->child_NE = NULL;
			else if(where == 2)
				cur->child_NW = NULL;
			else if(where == 3)
				cur->child_SW = NULL;
			else
				cur->child_SE = NULL;
		}
	}

	root->parent = NULL;
	root->child_SE = NULL;
	root->child_SW = NULL;
	root->child_NE = NULL;
	root->child_NW = NULL;
	root->starIndex = -1;
	root->starNumLower = 0;
	root->massPosX = 0;
	root->massPosY = 0;
	root->nodeWidth = boundWidth;
	root->centerX = boundMinX + boundWidth/2;
	root->centerY = boundMinY + boundWidth/2;
}

void applyForce(int indexStart, int indexEnd, long id) {

	//pthread_mutex_lock(&mutex);
	stars[indexStart].x = tempstars[indexStart].x;
	stars[indexStart].y = tempstars[indexStart].y;
	boundMinXArr[id] = stars[indexStart].x;
	boundMinYArr[id] = stars[indexStart].y;
	boundMaxXArr[id] = stars[indexStart].x;
	boundMaxYArr[id] = stars[indexStart].y; 
	for(int starIndex = indexStart+1; starIndex < indexEnd; starIndex++) {
		stars[starIndex].x = tempstars[starIndex].x;
		stars[starIndex].y = tempstars[starIndex].y;
		if(stars[starIndex].x < boundMinXArr[id])
			boundMinXArr[id] = stars[starIndex].x;
		if(boundMaxXArr[id] < stars[starIndex].x)
			boundMaxXArr[id] = stars[starIndex].x;
		if(stars[starIndex].y < boundMinYArr[id])
			boundMinYArr[id] = stars[starIndex].y;
		if(boundMaxYArr[id] < stars[starIndex].y)
			boundMaxYArr[id] = stars[starIndex].y;
	}
	//pthread_mutex_unlock(&mutex);
}

void computeForce(int indexStart, int indexEnd) {
	double Fx = 0;
	double Fy = 0;

	//pthread_mutex_lock(&mutex);
	for(int starIndex = indexStart; starIndex < indexEnd; starIndex++) {
		Fx = 0;
		Fy = 0;
		
		traverse(starIndex, Fx, Fy, 0, root, 0); //update force, argument is 0
		stars[starIndex].vx = stars[starIndex].vx + Fx;
		tempstars[starIndex].x = stars[starIndex].x + stars[starIndex].vx*timeOfStep;
		stars[starIndex].vy = stars[starIndex].vy + Fy;
		tempstars[starIndex].y = stars[starIndex].y + stars[starIndex].vy*timeOfStep;
		
	}
	//pthread_mutex_unlock(&mutex);
}

void traverse(int starIndex, double &Fx, double &Fy, int forceORpos, bhTree *now, int level) {
	double disX, disY, d, s;

	if(forceORpos == 0 && now != NULL) {
		disX = now->massPosX - stars[starIndex].x;
		disY = now->massPosY - stars[starIndex].y;
		d = sqrt(disX*disX + disY*disY);
		s = now->nodeWidth;
	}

	if((forceORpos == 0) && (now->starIndex == -1) && (s/d < theta)) {
		double r_3 = sqrt(disX*disX + disY*disY)*(disX*disX + disY*disY) + E;
		assert (abs(r_3) > 0.0);
		Fx += GMMS*(now->starNumLower)*disX/r_3;
		Fy += GMMS*(now->starNumLower)*disY/r_3;
		return;
	}

	if((forceORpos == 0) && (now->starIndex != -1)) {
		if(now->starIndex != starIndex) {
			double r_3 = sqrt(disX*disX + disY*disY)*(disX*disX + disY*disY) + E;
			assert (abs(r_3) > 0.0);
			Fx += GMMS*disX/r_3;
			Fy += GMMS*disY/r_3;
		}
		return;
	}

	if(forceORpos == 1) {
		double halfWidth = now->nodeWidth/2;
		double swX = now->centerX - halfWidth;
		double swY = now->centerY - halfWidth;
		double seX = now->centerX + halfWidth;
		double seY = now->centerY - halfWidth;
		double nwX = now->centerX - halfWidth;
		double nwY = now->centerY + halfWidth;
		double neX = now->centerX + halfWidth;
		double neY = now->centerY + halfWidth;
		double upX = now->centerX;
		double upY = now->centerY + halfWidth;;
		double downX = now->centerX;
		double downY = now->centerY - halfWidth;
		double leftX = now->centerX - halfWidth;
		double leftY = now->centerY;
		double rightX = now->centerX + halfWidth;
		double rightY = now->centerY;
		if(!strcmp(xwindowOpen, "enable")) {
			IOStart = gettime();
			XDrawLine(display, buffer, gc, (swX-xmin)*cor, (swY-ymin)*cor, (nwX-xmin)*cor, (nwY-ymin)*cor);
			XDrawLine(display, buffer, gc, (swX-xmin)*cor, (swY-ymin)*cor, (seX-xmin)*cor, (seY-ymin)*cor);
			XDrawLine(display, buffer, gc, (neX-xmin)*cor, (neY-ymin)*cor, (nwX-xmin)*cor, (nwY-ymin)*cor);
			XDrawLine(display, buffer, gc, (neX-xmin)*cor, (neY-ymin)*cor, (seX-xmin)*cor, (seY-ymin)*cor);
			if(now->starIndex == -1) {
				XDrawLine(display, buffer, gc, (upX-xmin)*cor, (upY-ymin)*cor, (downX-xmin)*cor, (downY-ymin)*cor);
				XDrawLine(display, buffer, gc, (leftX-xmin)*cor, (leftY-ymin)*cor, (rightX-xmin)*cor, (rightY-ymin)*cor);
			}
			IOTime = IOTime + gettime() - IOStart;
		}
	}

	if(forceORpos == 1 && now->starIndex == -1) {
		now->massPosX/=now->starNumLower;
		now->massPosY/=now->starNumLower;
	}

	if(forceORpos == 1 && now->starIndex != -1) 	
		return;

	if(now->child_SW != NULL) 
		traverse(starIndex, Fx, Fy, forceORpos, now->child_SW, level+1);

	if(now->child_SE != NULL) 
		traverse(starIndex, Fx, Fy, forceORpos, now->child_SE, level+1);

	if(now->child_NW != NULL) 
		traverse(starIndex, Fx, Fy, forceORpos, now->child_NW, level+1);

	if(now->child_NE != NULL)
		traverse(starIndex, Fx, Fy, forceORpos, now->child_NE, level+1);

	return;
}

void buildTree(int indexStart, int indexEnd) {
	buildTreeStart = gettime();
	int fall;
	bhTree *now;
	bhTree *nowParent;

	for(int starIndex = indexStart; starIndex < indexEnd; starIndex++) {
		now = root;
		nowParent = NULL;
		while(now != NULL) {
			//printf("whereToFall - starIndex %d assignX: %lf massX: %lf\n",starIndex, stars[starIndex].x, now->massPosX);
			fall = whereToFall(starIndex, &now, &nowParent);
			if((now != NULL && now->starIndex != -1) || now == NULL)
				break;
		}
		if(now != NULL && now->starIndex != -1) { //此層有star，要一起往下掉
			int oriStarIndex = now->starIndex;
			now->massPosX = 0;
			now->massPosY = 0;
			int starFall = whereToFall(starIndex, &now, &nowParent);
			now = nowParent;
			nowParent = now->parent;
			int oriStarFall = whereToFall(oriStarIndex, &now, &nowParent);
			while(starFall == oriStarFall) {
				create_node(starFall, -1, &now, &nowParent);
				starFall = whereToFall(starIndex, &now, &nowParent);
				now = nowParent;
				nowParent = now->parent;
				oriStarFall = whereToFall(oriStarIndex, &now, &nowParent);
			}
			create_node(starFall, starIndex, &now, &nowParent);
			create_node(oriStarFall, oriStarIndex, &now, &nowParent);

		} else { //此層目前為NULL，沒有star，可以建造node
			create_node(fall, starIndex, &now, &nowParent);
		}
	}
	buildTreeTime = buildTreeTime + gettime() - buildTreeStart;
}

void create_node(int fall, int starIndex, bhTree **now, bhTree **nowParent) {

	*now = (bhTree *)malloc(sizeof(bhTree));
	(*now)->parent = *nowParent;
	(*now)->child_SE = NULL;
	(*now)->child_SW = NULL;
	(*now)->child_NE = NULL;
	(*now)->child_NW = NULL;
	(*now)->nodeWidth = (*nowParent)->nodeWidth/2;
	(*now)->starIndex = starIndex;
	(*now)->starNumLower = 0;
	if(starIndex != -1) {
		(*now)->massPosX = stars[starIndex].x;
		(*now)->massPosY = stars[starIndex].y;
	} else {
		(*now)->massPosX = 0;
		(*now)->massPosY = 0;
	}
	if(fall == 1) {
		(*now)->centerX = (*nowParent)->centerX + (*now)->nodeWidth/2;
		(*now)->centerY = (*nowParent)->centerY + (*now)->nodeWidth/2;
		(*nowParent)->child_NE = *now;
	}
	else if(fall == 2) {
		(*now)->centerX = (*nowParent)->centerX - (*now)->nodeWidth/2;
		(*now)->centerY = (*nowParent)->centerY + (*now)->nodeWidth/2;
		(*nowParent)->child_NW = *now;
	}
	else if(fall == 3) {
		(*now)->centerX = (*nowParent)->centerX - (*now)->nodeWidth/2;
		(*now)->centerY = (*nowParent)->centerY - (*now)->nodeWidth/2;
		(*nowParent)->child_SW = *now;
	}
	else {
		(*now)->centerX = (*nowParent)->centerX + (*now)->nodeWidth/2;
		(*now)->centerY = (*nowParent)->centerY - (*now)->nodeWidth/2;
		(*nowParent)->child_SE = *now;
	}
}

int whereToFall(int starIndex, bhTree **now, bhTree **nowParent) {
	int fall = 0;
	double assignX = stars[starIndex].x;
	double assignY = stars[starIndex].y;

	(*now)->starIndex = -1;
	//printf("massX: %lf\n", (*now)->massPosX);
	(*now)->massPosX+=assignX;
	(*now)->massPosY+=assignY;
	//printf("After - massX: %lf\n", (*now)->massPosX);
	(*now)->starNumLower++;
	//printf("starIndex: %d cur->starNumLower: %d\n", starIndex, cur->starNumLower);
	*nowParent = *now;

	/*決定要往哪邊走，分為四個象限*/
	if(assignX < (*now)->centerX) {
		if(assignY < (*now)->centerY) { 
			fall = 3;
			*now = (*now)->child_SW;
		}
		else { 
			fall = 2;
			*now = (*now)->child_NW;
		}
	}
	else {
		if(assignY < (*now)->centerY) {
			fall = 4;
			*now = (*now)->child_SE;
		}
		else {
			fall = 1;
			*now = (*now)->child_NE;
		}
	}

	return fall;
}

void initGraph(int width, int height) {
	/* open connection with the server */ 
	display = XOpenDisplay(NULL);
	if(display == NULL) {
		fprintf(stderr, "cannot open display\n");
		exit(1);
	}

	screen = DefaultScreen(display);

	/* set window position */
	int x = 0;
	int y = 0;

	/* border width in pixels */
	int border_width = 0;

	/* create window */
	window = XCreateSimpleWindow(display, RootWindow(display, screen), x, y, width, height, border_width, BlackPixel(display, screen), WhitePixel(display, screen));
	
	/* create graph */
	XGCValues values;
	long valuemask = 0;
	
	gc = XCreateGC(display, window, valuemask, &values);
	//XSetBackground (display, gc, WhitePixel (display, screen));
	XSetForeground (display, gc, BlackPixel (display, screen));
	XSetBackground(display, gc, 0X0000FF00);
	XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);
	
	/* map(show) the window */
	XMapWindow(display, window);
	XSync(display, 0);

	buffer = XCreatePixmap(display, window, xwindowLength, xwindowLength, 24);

	/* draw rectangle */
	XSetForeground(display,gc,BlackPixel(display,screen));
	XFillRectangle(display,window,gc,0,0,width,height);
	XFlush(display);
}

void draw(int x, int y) {
	/* draw point */
	XSetForeground(display,gc,WhitePixel(display,screen));
	XDrawPoint (display, buffer, gc, x, y);
}
