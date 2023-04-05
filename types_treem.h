/* defs.h */

//#ifdef EXCESS_TYPE_LONG
///typedef unsigned long excessType;
//#else
typedef unsigned long long int excessType; /* change to double if not supported */
//#endif

typedef unsigned long cType;
typedef unsigned int sType;

typedef 
   struct edgeProp
{

   cType endNode;
   cType cap; //边的值
   cType w; /*weight*/
   cType avgCV; /*average CV of cut set where this edge belongs to */
   long tmp;

   struct edgeProp* rev; //reverse 边反过来

}edgeP;


typedef  /* node */
   struct nodeProp
{
   edgeP* edges;
   cType maxEdges;
   cType nIdx; //节点的边的个数
   cType totalCap;
   
   sType* orderedEdges; //给每个边一个序号

} nodeP;

typedef
struct NodePropExtra_
{
   cType fa;
   cType dep;
   long cv;
   short s;


} NodePropExtra;


typedef
struct NodePropArr_{

   cType * pfa;
   cType * pdep;
   long * pcv;
   long * poh; //when traversing: record the edges' capacity connected to f by tree of a child n 遍历时：通过子n的树记录连接到f的边的容量
               //when traversed: record the new updated cv 遍历时：记录新更新的cv
   long * pcof; //record the remove of n, bring to the change of father' cv 记录n的去除，带来父亲简历的变化
   short * ps;
   cType * pacc_upid; // up node id 上行节点id
   long* pacc_upmincv; // min of [currnet node, upnode) [currnet节点，upnode的最小值
   cType * pacc_pos_upmincv; // the diff of depth of nearest min value 最近最小值的深度差
   cType * pacc_jointid; //the nearest joint id; 最近的关节id
   long* pacc_jointmcv; //the min of [cur, joint node)  [cur，joint node）的最小值
   
} NodePropArr;

