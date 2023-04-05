/* Maximum flow - highest lavel push-relabel algorithm */
/* COPYRIGHT C 1995, 2000 by IG Systems, Inc., igsys@eclipse.net */
  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <values.h>
#include <math.h>

#include "types_treem.h"  /* type definitions */
#include "parser_treem.c" /* parser */
#include "timer.c"        /* timing routine */

#define MAX_LONG LONG_MAX

#define min(s, t) (s) < (t) ? (s) : (t)
#define max(s, t) (s) > (t) ? (s) : (t)

/////////////////////////////////////////
///////////////////////////////////////// The definition of data structure
/////////////////////////////////////////

//data structure that holds graph data
typedef
struct GraphData_{

  long N, M;
  nodeP *nodes; //节点数组,每个节点都是一个节点及其所有边到ID比它大的节点信息

} GraphData;


//data structure that holds data for randomrization
typedef
struct RandomData_{
  cType *randNums;
  cType randNumIdx;
} RandomData;


//data structure that holds preprocessing data
typedef
struct PreprocData_{

  //graph data
  GraphData *gd;

  //holds data in all passes
  NodePropArr* allResults;

  //pre-generated data  预生成数据
  RandomData* rd;

  //<BEING> hold the hot data in current pass
  cType *gpfa;    //[i]是节点i的父节点id
  cType *gpdep;  //[i]是节点i的遍历深度
  long *gpcv;  //[i]是节点i的CV
  short *gps; //[i]是节点的状态
  cType *gpaccup; //[i]是算法3中节点i的sh
  long *gpaccmcv;//[i]是算法3中节点i的mcv
  cType *gpaccposmcv;//对算法3的一点改进
  cType *gpaccjointid;//未用
  long *gpaccjointmcv;//未用
  //<END>

  cType *roots; //records root in each pass 记录每次通过的根

  int mode; //the traversing mode, as explained in the paper

  int P; // P% percent of total passes in mode 1, the remaining in mode 2
  int total; //number of total passes 总通过次数
  int SPAN_LEN; //the length of a ancestor span, used for acceleration in traversal trees 祖先跨度的长度,用于遍历树中的加速

} PreprocData;




///////////////////////////////////////////
///////////////////////////////////////////The definition of functions
///////////////////////////////////////////


//The function for allocation
void *walloc(unsigned int num, unsigned int size)
{
  void *ptr = calloc(num, size);
  assert(ptr != NULL);
  return ptr;
}

//The function for randomization
cType mrand(RandomData* rd)
{
  return rd->randNums[rd->randNumIdx++];
}

RandomData* initrand(cType len)
{
  RandomData *rd = walloc(1, sizeof(RandomData));

  srand((int)(timer() * 1000));

  rd->randNums = (cType *)walloc(len, sizeof(cType));
  rd->randNumIdx = 0;

  for (int i = 0; i < len; i++)
  {
    rd->randNums[i] = (cType)rand();
  }

  return rd;
}

/////////////////////////The two function for heap sorting edges  堆排序的两个函数,在算法2中用于对节点边缘排序
void HeapAdjustDown(sType *idx, edgeP * edges ,int start,int end)  
{  
    sType tempIdx = idx[start];  
    //printf("start是%d\n",tempIdx);
    //printf("当前边的temp:%d\n",edges[tempIdx].tmp);
    int i = 2*start+1;      
    
    // assert(idx[0] != idx[3]);
    
    while(i<=end)  
    {  
        if(i+1<=end && edges[idx[i+1]].tmp > edges[idx[i]].tmp )    
            i++;  

        if(edges[idx[i]].tmp <= edges[tempIdx].tmp )   
            break;  

        idx[start] = idx[i];

        start = i;  
        i = 2*start+1;  
    }  

    idx[start] = tempIdx;  
}  
  
void HeapSort(sType *idx, edgeP * edges, int len)  
{  
  //printf("len是%d \n",len);
    int i;  
    for(i=(len-1)/2;i>=0;i--){  
        //printf("i是%d \n",i);
        HeapAdjustDown(idx,edges,i,len-1);  
    }

    for(i=len-1;i>0;i--)
    {  
        // printf("swap 0 with %d \n",i);
        sType temp = idx[i];  
        idx[i] = idx[0];  
        idx[0] = temp;  

        HeapAdjustDown(idx,edges,0,i-1);  
    }  

}  
  
/////////////////////The function to sort edges using capacity 使用容量对边进行排序
void deOrderEdgeByRandomCap(nodeP *np,PreprocData *pd)
{

  cType cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    pedges[i].tmp = -1*mrand(pd->rd) % pedges[i].cap;
    //printf("pedges[%d].tmp:%d \n",i,pedges[i].tmp);
  }

    assert(cnt<4 || idxs[2]!=idxs[3]);
    HeapSort(idxs,pedges,cnt);
    assert(cnt<4 || idxs[2]!=idxs[3]);  
}

///////////////////The function to sort edges using the value of currently minimal cut one edge belongs to  使用当前最小割的值对边排序
void aOrderEdgeByAvgCV(nodeP *np,PreprocData *pd)
{
  if (np->nIdx == 0)
  {
    return;
  }
  int cnt = np->nIdx;

  sType *idxs = np->orderedEdges;
  edgeP *pedges = np->edges;

  for (int i = 0; i < cnt; i++)
  {
    edgeP *pp = pedges +i;
    long acv = pp->avgCV;
    if(acv == 0 ){
      pp->tmp = MAX_LONG;
    }
    else{
      pedges[i].tmp = mrand(pd->rd) % acv; 
    }
  }

    HeapSort(idxs,pedges,cnt);
}


/*
  the function to add more data to a traversal tree to accelerate the searching in the tree
  The idea is to precalcuate minimal cv value of a span of nodes in a traversal tree, e.g., when SPAN_LEN = 100, and a node has depth of 200, then the algorithm will pre-calculate the minimal cv value of the nodes between the node (dep=200) and an acestor(dep=101)
  upid is the id of the last SPAN node, mcv is the min cv among all previous nodes in the recent SPAN
  lastDepMCV is the depth of the node depth that has the minimal cv in the span
  lastJointNodeId is the last ancestor node id that has more than one child nodes
  lastJointMCV is the cv of lastJoineNodeId
向遍历树中添加更多数据以加快树中搜索的功能
其思想是预先计算遍历树中节点跨度的最小cv值，例如，当span_LEN=100，并且节点的深度为200时，
算法将预先计算节点（dep=200）和acestor（dep=101）之间节点的最小cw值
upid是最后一个SPAN节点的id，mcv是最近SPAN中所有先前节点的最小cv
lastDepMCV是跨度中具有最小cv的节点深度的深度
lastJointNodeId是具有多个子节点的最后一个祖先节点id
lastJointMCV是lastJoineNodeId的cv
*/
void buildAcc(PreprocData *pd, cType curN, cType upid, long mcv, cType lastDepMCV,cType lastJointNodeId, long lastJointMCV)
{
  nodeP* nodes = pd->gd->nodes;
  assert(curN <= pd->gd->N && curN >= 1);
  int cnt = (nodes + curN)->nIdx;
  //printf("当前节点为%d cnt is %d\n",curN,cnt);
  edgeP *pedges = (nodes + curN)->edges;

  long curCV = pd->gpcv[curN];
  //printf("curCV is %d mcv is %d\n",curCV,mcv);
  cType curDep = pd->gpdep[curN];
  //printf("curDep is %d\n",curDep);

  assert(curDep == 0 || curCV >0);
//logic for span
  mcv = min(mcv, curCV);//一个span内最小cv值
  if(mcv == curCV){
    lastDepMCV = curDep;
  }
  //printf("lastDepMCV is %d\n",lastDepMCV);

  pd->gpaccmcv[curN] = mcv;
    //printf("pd->gpaccmcv[curN] is %d\n",pd->gpaccmcv[curN]);
  pd->gpaccup[curN] = upid;
    //printf(" 当前span的跟节点pd->gpaccup[curN] is %d\n", pd->gpaccup[curN]);
  pd->gpaccposmcv[curN] = curDep - lastDepMCV;
  //printf("与最小mv节点深度差pd->gpaccposmcv[curN] is %d\n",pd->gpaccposmcv[curN]);

  if (curDep % pd->SPAN_LEN == 0)
  {
    upid = curN;
    mcv = MAX_LONG;
    lastDepMCV = 0; //doesn't matter, will be udpated in the subsequent call
  }

//-logic for joint node
  cType childCnt = pd->gpaccjointid[curN];
  //printf("childCnt is %d\n",childCnt);
  lastJointMCV = min(lastJointMCV,curCV);
  //printf("全局最小mv值lastJointMCV is %d\n",lastJointMCV);
  pd->gpaccjointmcv[curN] = lastJointMCV; 
  pd->gpaccjointid[curN] = lastJointNodeId;
  //printf("总跟节点pd->gpaccjointid[curN] is %d\n\n",pd->gpaccjointid[curN]);

  assert(pd->gpdep[curN] == 0 || pd->gpaccjointmcv[curN] > 0.1);
  
  //gpaccjointid[curN] is updated by markcut() to contain the number of traversing children
  if(childCnt > 1 ){
    lastJointNodeId = curN;
    lastJointMCV = MAX_LONG;
  }
  else{
    //do nothing
  }


  while (cnt > 0)
  {
    if (pd->gpfa[pedges->endNode] == curN)
    {
      buildAcc(pd,pedges->endNode, upid, mcv,lastDepMCV,lastJointNodeId,lastJointMCV);
    }
    pedges++;
    cnt--;
  }

  // assert(childCnt == 0);
}

//////////////////////////////////The function to traverse the graph for one pass, i.e. the checkNode function of Algorithm 2 in the paper //遍历一次图像,也就是文中checkNode函数
void markCut(cType curN, PreprocData *pd)
{
  nodeP* nodes = pd->gd->nodes;
  //printf("nodes  maxEdges:%d nIdx:%d totalCap:%d orderedEdges:%d",nodes->maxEdges,nodes->nIdx,nodes->orderedEdges,nodes->totalCap);
  //printf("curN is %d ",curN);

  assert(curN <= pd->gd->N && curN >= 1);

  assert((nodes + curN)->nIdx > 0);

  short *curS = pd->gps + curN;
 //printf("pd->gps is %d",pd->gps);
//printf("curS is %d ",pd->gps + curN);
  assert(*curS == 0);

  long *curCV = pd->gpcv + curN;
  //printf("curCV is %d\n",pd->gpcv + curN);
  cType *curDep = pd->gpdep + curN;

  *curS = 1;
  *curCV = 0;
  nodeP *np = nodes + curN;
  edgeP *pedges = np->edges;
  int cnt = np->nIdx;
  //printf("cnt is %d\n",cnt);
  //printf("初始np  maxEdges:%d nIdx:%d totalCap:%d orderedEdges:%d\n",np->maxEdges,np->nIdx,np->orderedEdges,np->totalCap);
 // printf("初始pedges endNode:%d cap:%d w:%d  avgCV:%d tmp:%d\n",np->edges->endNode,np->edges->cap,np->edges->w,np->edges->avgCV,np->edges->tmp);

  if (pedges == NULL)
  {
    *curS = 2;
    return;
  }

  if (np->orderedEdges == NULL)
  {
    np->orderedEdges = (sType *)walloc(cnt + 1, sizeof(sType));
    for (int i = 0; i < cnt; i++)
    {
      np->orderedEdges[i] = i;
      //printf("np orderedEdges[%d]:%d\n",i,np->orderedEdges[i]);
    }
  }
  //printf("np orderedEdges:%d\n",np->orderedEdges);

  long cap;
  sType *idxs = np->orderedEdges;

  if (pd->mode == 1)
  {
    deOrderEdgeByRandomCap(np,pd); //
  }
  else if (pd->mode == 2)
  {
    aOrderEdgeByAvgCV(np,pd);
  }
  //printf("排序后np  maxEdges:%d nIdx:%d totalCap:%d orderedEdges:%d\n",np->maxEdges,np->nIdx,np->orderedEdges,np->totalCap);
  //printf("排序后pedges endNode:%d cap:%d w:%d  avgCV:%d tmp:%d\n\n",np->edges->endNode,np->edges->cap,np->edges->w,np->edges->avgCV,np->edges->tmp);

  for (int ni = 0; ni < cnt; ni++)
  {
    // nodeP* znp = nodes+eh->endNode;

    edgeP *eh = pedges + idxs[ni];
    cType zn = eh->endNode;
   //printf("eh pedges endNode:%d cap:%d w:%d  avgCV:%d tmp:%d\n",eh->endNode,eh->cap,eh->w,eh->avgCV,eh->tmp);
    assert(zn != 0);
    assert(zn != curN);

    short zs = pd->gps[zn]; //当前节点状态
    
   //printf("zs:%d \n",zs);

    assert(!(pd->gpfa[zn] == curN && zs == 2));


    if (zs == 1) //说明该节点以遍历
    {
      cap = eh->cap;
      *curCV += cap;
      pd->gpcv[zn] -= cap;
      //printf("当前%d节点cv:%d \n",zn,pd->gpcv[zn]);

    }
    else if (zs == 0) //递归该节点
    {
      pd->gpfa[zn] = curN;
      pd->gpdep[zn] = *curDep + 1;
      pd->gpaccjointid[curN] ++;
      markCut(zn,pd);
      assert(pd->gpdep[zn] == pd->gpdep[curN] + 1);
      *curCV += pd->gpcv[zn];
        //printf("当前%d节点cv:%d \n",zn,pd->gpcv[zn]);

    }
    else
    {
      //bypass, no need to handle
    }
  }
//printf("没有未访问节点\n");

  if(pd->mode == 1){
//update ver 2 w,according to
    for (int ni = 0; ni < cnt; ni++)
    {
      // nodeP* znp = nodes+eh->endNode;

      edgeP *eh = pedges + idxs[ni];
      cType zn = eh->endNode;
      // nodeP *znp = nodes+zn;
    //printf("pedges endNode:%d cap:%d w:%d  avgCV:%d tmp:%d\n",eh->endNode,eh->cap,eh->w,eh->avgCV,eh->tmp);

      assert(zn != 0);
      assert(zn != curN);
      short zs = pd->gps[zn];
    

      //progate weight to curN's edges 将权重累加到curN的边
      if (zs == 1 && pd->gpdep[zn] != *curDep - 1)
      {
          cType weight = eh-> w;
          if(eh->avgCV == 0){
            eh->avgCV = MAX_LONG;
          }
          //printf("eh->avgCV:%d  *curCV:%d\n",eh->avgCV, *curCV);
          eh->avgCV = min(eh->avgCV, *curCV);//((eh->avgCV) * weight + *curCV)/(weight+1);
          eh->w = weight+1;

          edgeP *reh = eh->rev;

          if(reh->avgCV == 0){
            reh->avgCV = MAX_LONG;
          }

          weight = reh-> w;
          reh->avgCV = min(reh->avgCV, *curCV);//((reh->avgCV) * weight + *curCV)/(weight+1);
          reh->w = weight+1;

      }
      
    }
  }

    
  assert(pd->gpdep[curN] == 0 || *curCV >0);

  *curS = 2;

}

////////////////////////////////The function to obtain min-cut value of given node pair, i.e., Algorithm 3 in the paper 获取给定节点对的最小割
long solveMaxFlowAccVER4(long minCandi, NodePropArr np, cType s, cType t, int SPAN_LEN)
{
  cType *pDep = np.pdep;
  //printf("pDep is %d\n",*pDep);
  long *pCV = np.pcv;
  //printf("pCV is %d\n",*pCV);
  cType *pFa = np.pfa;
  //printf("pFa is %d\n",*pFa);
  cType *paccup = np.pacc_upid;
  //printf("上行节点id paccup is %d\n",*paccup);
  long *paccmcv = np.pacc_upmincv;
  //printf("paccmcv is %d\n",*paccmcv);
  cType *paccposmcv = np.pacc_pos_upmincv;
  //printf("paccposmcv is %d\n",paccposmcv);
  cType *jup = np.pacc_jointid;
  //printf("最近关节点jup is %d\n",*jup);
  long *jmcv = np.pacc_jointmcv;
  //printf("jmcv is %d\n\n",*jmcv);

  assert(s != t);
  if (pDep[s] < pDep[t])
  {
    cType tmp = s;
    s = t;
    t = tmp;
  }

  cType depT = pDep[t];
  //printf("pDeps is %d depT is %d\n",pDep[s],depT);

  assert(pDep[s] >= depT);

  long mcv = MAX_LONG;
  cType ups = paccup[s]; //找到上一个段尾
 //printf("上个段尾paccup is %d\n",ups);

  while (pDep[ups] > depT)  //如果pDep[ups] =depT就不会继续,pDep[ups] > depT此时ups的cv'肯定纳入计算
  {
    assert(pDep[ups] % SPAN_LEN == 0);
    mcv = min(mcv, paccmcv[s]);
    // printf("当前mcv是 %d s是%d\n",mcv,s);
    s = ups;
    ups = paccup[ups];
  } 
//   assert(mcv >100.0);
  assert(pDep[ups] <=depT && pDep[s] >= depT);

  cType upt = t;

  if (pDep[t] % SPAN_LEN != 0)
  {
    upt = paccup[t];
  }

  //assert(pDep[ups] == pDep[upt]);

  while (ups != upt)
  {
    mcv = min(mcv, min(paccmcv[s], paccmcv[t]));
    //printf("当前mcv is %d\n",mcv);
    s = ups;
    ups = paccup[ups];

    t = upt;
    upt = paccup[upt];

    assert(pDep[s] % SPAN_LEN == 0);
    //assert(pDep[t] == pDep[s]);
  }        
  
  assert(ups == upt);
  if(s == t){
    return mcv;
  }
  
  cType min_bound2 = min(paccmcv[s],paccmcv[t]);

  if(min_bound2 >= mcv || min_bound2 >= minCandi){
    //no need to search
    return mcv;
  }

  if(min_bound2 == paccmcv[s] && pDep[t]+paccposmcv[s] < pDep[s]){
    return mcv;
  }


  ///////////////////////check inside one SPAN
  //(1) we need to check whether s and t in the same line, i.e., t is the ancestor of s  检查s和t是否在同一行,即t是s的祖先
  //the only way is to apprach the depth of t and check 唯一办法是接近t的深度并检查

  if(pDep[s] != pDep[t]){

    if(pDep[s] < pDep[t]){
      cType tmp = s;
      s = t;
      t = tmp;
    }  

    cType jups = jup[s];
    cType depT = pDep[t];
    while(pDep[jups] > depT){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }

      s = jups;
      jups = jup[s];

    }

    if(pDep[jups] == depT){
      if(jups == t){
        return min(mcv,jmcv[s]);
      }
      else{
        //s and t in two lines
        assert(pDep[jups] == depT && pDep[s] > depT);
        goto STEP_CHECK_IN_TWOLINES;
      }
    }
    else{
      assert(pDep[jups] < depT && pDep[s] > depT);
      while(pDep[s] > depT){
        mcv = min(mcv,pCV[s]);
        if(mcv == min_bound2){
          return mcv;
        }
        s = pFa[s];        
      }

      if(s == t){
        return mcv;
      }

      assert(s!=t && pDep[s] == depT);
      //s and t in two lines.
      goto STEP_CHECK_IN_TWOLINES;

    }
  }
  else{
      //pDep[s] == pDep[t];
      if(s == t){
        return mcv;
      }    
      goto STEP_CHECK_IN_TWOLINES;
  }

STEP_CHECK_IN_TWOLINES:

  while (s != t)
  { 
    // assert(jmcv[s] > 0.1);
    // assert(jmcv[t] > 0.1);
    if(pDep[s] > pDep[t]){
      mcv = min(mcv, jmcv[s]);
      if(mcv == min_bound2){
        return mcv;
      }      
      s = jup[s];
    }
    else if(pDep[s] < pDep[t]){
      mcv = min(mcv, jmcv[t]);
      if(mcv == min_bound2){
        return mcv;
      }      
      t = jup[t];
    }
    else{
      mcv = min(mcv, min(jmcv[s],jmcv[t]));
      if(mcv == min_bound2){
        return mcv;
      }         
      s = jup[s];
      t = jup[t];
    }

  }

  return mcv;

}

//////////////////////function to load graph data
void loadGraphData(PreprocData *pd){
  pd->gd = walloc(1,sizeof(GraphData));  
  printf("c\nc hi_treem version 0.9\n");
  printf("c Copyright C by nsyncw, nsyncw@gmail.com\nc\n");

  parse(&(pd->gd->N), &(pd->gd->M), &(pd->gd->nodes));

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", pd->gd->N, pd->gd->M);
}

///////////////////function to initialize data structure
void initPreprocData(PreprocData *pd){
  pd->rd = NULL;
  pd->SPAN_LEN = (int)(sqrt(pd->gd->N));

  pd->roots = (cType *)walloc(pd->total + 2, sizeof(cType));
  pd->allResults = walloc(pd->total+2, sizeof(NodePropArr));

  NodePropArr * allResults = pd->allResults;
  cType len = pd->gd->N + 2;
  for (int i = 0; i < pd->total; i++)
  {
    allResults[i].pfa = (cType *)walloc(len, sizeof(cType));
    allResults[i].pdep = (cType *)walloc(len, sizeof(cType));
    allResults[i].pcv = (long *)walloc(len, sizeof(long));
    allResults[i].ps = (short *)walloc(len, sizeof(short));
    allResults[i].pacc_upid = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_upmincv = (long *)walloc(len, sizeof(long));
    allResults[i].pacc_pos_upmincv = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_jointid = (cType *)walloc(len, sizeof(cType));
    allResults[i].pacc_jointmcv = (long *)walloc(len, sizeof(long));

    memset(allResults[i].pfa, 0, len * sizeof(cType));
    memset(allResults[i].pdep, 0, len * sizeof(cType));
    memset(allResults[i].pcv, 0, len * sizeof(long));
    memset(allResults[i].ps, 0, len * sizeof(short));
    memset(allResults[i].pacc_upid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_upmincv, 0, len * sizeof(long));
    memset(allResults[i].pacc_pos_upmincv, 0, len * sizeof(cType));
    memset(allResults[i].pacc_jointid, 0, len * sizeof(cType));
    memset(allResults[i].pacc_jointmcv, 0, len * sizeof(long));
  }  
}


/////////////////////function to traverse the graph data for multiple times, i.e., Algorithm 2 in the paper
void preProc(PreprocData *pd){
  double tm;
  double totalProcTime = 0;
  NodePropArr *allResults = pd->allResults;
  //calculate total cap of one node
  cType root;
  for (int ipass = 0; ipass < pd->total; ipass++)
  {
    if(pd->rd != NULL){
      free(pd->rd);
      pd->rd = NULL;
    }
    pd->rd = initrand(pd->gd->M*2);    
    // printf("the %d times\n",i);
    pd->gpfa = allResults[ipass].pfa;
    pd->gpdep = allResults[ipass].pdep;
    pd->gpcv = allResults[ipass].pcv;
    pd->gps = allResults[ipass].ps;
    pd->gpaccup = allResults[ipass].pacc_upid;
    pd->gpaccmcv = allResults[ipass].pacc_upmincv;
    pd->gpaccposmcv = allResults[ipass].pacc_pos_upmincv;
    pd->gpaccjointid = allResults[ipass].pacc_jointid;
    pd->gpaccjointmcv = allResults[ipass].pacc_jointmcv;

    pd->mode = ipass < pd->P * pd->total / 100 ? 1 : 2;
    
    root = 1 + ((mrand(pd->rd) * mrand(pd->rd)) % pd->gd->N);

    pd->roots[ipass] = root;
    pd->gpdep[root] = 0;
    //printf("pass %d, randidx %ld, root is %ld\n",ipass, randNumIdx,root);
    // printf("pass %d, befor marcut root is %ld\n",ipass,root);
    // printf("root fa %ld\n",allResults[ipass].pfa[root]);
    tm = timer();
    markCut(root,pd);

    pd->gpcv[root] = MAX_LONG;
     //printf("after marcut  root is%d pd->gpcv[root]:%d\n",root,pd->gpcv[root]);

    buildAcc(pd, root, root, MAX_LONG, pd->gpdep[root],root,MAX_LONG);

    totalProcTime += timer() - tm;
    printf("c proctime for onepass: %10.06f\n", timer() - tm);
    if (ipass % 10 == 0)
    {
      printf("c the %d passes\n", ipass);
    }
  }

  printf("c preprocess times %10.6f\n", totalProcTime);

}


//////////////////////////function to calculate multiple random node pairs, i.e., the calling of Algoirthm 3 for multiple times
void calcuRandomPairs(PreprocData *pd){
  int N = pd->gd->N;
  double totalTime = 0;
  long mv = MAX_LONG;

  double curTime = 0;
  cType ns, nt;
  int MAX_FLOW[N+1];
  int sum;
   curTime = timer();
  for(int i=1;i<=N;i++)
  {
    ns = i;
    sum = 0;
    for(int x=1;x<=N;x++)
    {
      nt = x;
      mv = MAX_LONG;
      if(ns == nt)
      {
        continue;
      }
      for (int j = 0; j < pd->total; j++)
      {
        long tmp = solveMaxFlowAccVER4(mv, pd->allResults[j], ns, nt,pd->SPAN_LEN);
        if (mv > tmp)
        {
          mv = tmp;
        }
        //printf("此次循环mv值为:%d,ns和nt分别为%d,%d\n",mv,ns,nt);
      }
      //printf("最终mv值为:%d,ns和nt分别为%d,%d\n",mv,ns,nt);
      sum = sum+mv;
    }
    MAX_FLOW[i] = sum;
  }
      curTime = timer() - curTime;
      totalTime += curTime;
      int ANF=0;
  for(int i=1;i<=N;i++)
  {
    printf("%d节点的总网络流为:%d\n",i,MAX_FLOW[i] );
    ANF=ANF+MAX_FLOW[i];
  }
  printf("ANF is %d\n", ANF/(N-1));
  //printf("c run ok! average time %10.6f\n", totalTime / N);
  //printf("c run ok! time %10.6f\n", totalTime );
  }


