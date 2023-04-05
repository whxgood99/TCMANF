#include "hi_treem.c"

int main(argc, argv)

    int argc;
char *argv[];

{
    PreprocData *pd = walloc(1,sizeof(PreprocData));
    pd->P = 30; // P% percent of total passes in mode 1, the remaining in mode 2
    pd->total =10000; //number of total passes
    double curTime = 0;
    double totalTime = 0;
    curTime = timer();

    loadGraphData(pd); //load graph data from standard input

    initPreprocData(pd); //init data structure

    preProc(pd); // preproc by traversing the graph for pd->total times，the Lines 1-6 in Algorithm 1 of the paper

    calcuRandomPairs(pd); // randomly choose 100 node pairs and calcu their min-cut and output. For each pair, it is the logics in Lines 7-10 in Algorithm 1.
    curTime = timer() - curTime;
    totalTime += curTime;
    printf("c run ok! time %10.6f\n", totalTime );
    exit(0);

}
