#pragma once

// only for eclipse/VS for syntax highlighting
#if ! defined(__OPENCL_VERSION__) &&  ! defined(__CUDACC__)
/* Define out keywords causing errors */
#define __kernel
#define __global
#define __constant
#define __local
#define __private
#define CLK_LOCAL_MEM_FENCE 32
#define align /*__declspec(align(8)) */
#else
#define align /*__attribute__ ((aligned))*/
#endif

#define chunk_size 256u
#define T 256u
#define G 43007.1

struct BBox
{
	union {
		FLOAT3 box[2];
		FLOAT boxArr[2][4];
	};
};

//#pragma pack(push, 4)
struct Particle
{
	union {
		FLOAT3 pos;
		struct{
			FLOAT posArr[3];
			FLOAT mass;
		};
	};

	union {
		FLOAT3 acc;	//force acting on the particle, calculated with tree walk
		struct{
			FLOAT accArr[3];
			FLOAT totAcc;
		};
	};
	
	union {
		FLOAT3 vel; //velocity & timestep of the particle
		struct{
			FLOAT velArr[3];
			//FLOAT dt; //not needed
			//UINT timebin; //timebin of the particle, used later for individual timesteps instead of dt
			UINT id;
		};
	};

//	UINT id; moved to separate array to shrink particle struct and solve alignment issues
};
//#pragma pack(pop)

//#pragma pack(push, 16)
typedef UINT NodeId;

struct Node
{
	//force tree
	union {
		FLOAT3 center_of_mass;		//center of mass
		struct{
			FLOAT cmArr[3];
			FLOAT mass;				//total mass in the node
		};	
	};
	

	union {
		FLOAT3 center_geometric;	//geometric center of the node's bounding box	
		struct{
			FLOAT cgArr[3];
			FLOAT l;					//extension
		};
	};

	//tree algorithm
	struct BBox bounding_box;	//bounding box of the node (min max for each dimension)

	UINT level; 		//level of the node in the tree
	UINT size; 			// size of tree (in nodes), including the subtree underneath it
	UINT address;		// index in final tree array

	UINT particlesLow;	// index of first particle belonging to that node
	UINT particlesHigh;	// index of last particle belonging to that node

//	UINT small;
//	std::vector<Chunk> chunks;

	NodeId left_child;
	NodeId right_child;

	UINT split_dim;		// dimension along node is split in child nodes. Accelerates sorting into childnodes
} align;

struct KdNode
{
	//force tree
	union {
		FLOAT3 center_of_mass;		//center of mass
		struct{
			FLOAT cmArr[3];
			FLOAT mass;				//total mass in the node
		};	
	};
	

	union {
		FLOAT3 center_geometric;	//geometric center of the node's bounding box	
		struct{
			FLOAT cgArr[3];
			FLOAT l;					//extension
		};
	};

	UINT size; 			// size of tree (in nodes), including the subtree underneath it

	UINT leaf;	// is leave node

//	UINT small;
//	std::vector<Chunk> chunks;

	NodeId left_child;
	NodeId right_child;
} align;


struct Chunk {
	UINT parentNode;
	UINT start;
	UINT end;
	UINT pad;
};

struct Tree
{
	struct Node* nodelist;
};
//#pragma pack(pop)
