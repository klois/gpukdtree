#ifndef KDTREE_H
#define KDTREE_H


#include <vector>
#include <iostream>
#include <cfloat>
#include <string.h>

#define FLOAT 		double
#define UINT		unsigned
#define FMAX		FLT_MAX

const UINT chunk_size = 256;
const UINT T = 256;
const FLOAT G = 43007.1;

class Node;
class Tree;

class BBox
{
public:
	FLOAT box[2][3];	//bounding box
	
	BBox() 
	{
		box[0][0] = box[0][1] = box[0][2] = FMAX;
		box[1][0] = box[1][1] = box[1][2] = 0.0;
	};
	//TODO: copy constr. ??

	FLOAT getVolume() const {
		return (box[1][0] - box[0][0]) * (box[1][1] - box[0][1]) * (box[1][2] - box[0][2]);
	}
	
	void getCenter(FLOAT center[3]) const
	{
		center[0] = (box[1][0] + box[0][0])/2.0;
		center[1] = (box[1][1] + box[0][1])/2.0;
		center[2] = (box[1][2] + box[0][2])/2.0;
	}
	
	friend std::ostream& operator<<(std::ostream &stream, BBox &bbox)
	{
		stream << "bbox: " << bbox.box[0][0] << " " << bbox.box[1][0] << " : "
			<< bbox.box[0][1] << " " << bbox.box[1][1] << " : "
			<< bbox.box[0][2] << " " << bbox.box[1][2] << std::endl;
		
		return stream;
	}

};

struct Particle
{
	UINT id;

	FLOAT pos[3];
	FLOAT mass;
	
	FLOAT acc[3];	//force acting on the particle, calculated with tree walk
	FLOAT totAcc;
	
	//time integration
	FLOAT dt; //timestep of this particle
	//UINT timebin; //timebin for this particle ...(just for individual timesteps, then we don't need dt anymore)
	
	FLOAT vel[3]; //velocity
	

	Particle() {
		acc[0] = acc[1] = acc[2] = 0;
		totAcc = 0;
	}
};

struct Chunk
{
	std::vector<Particle *> particles;
	BBox bounding_box;
	Node *node;					//pointer to the node the chunk belongs to
};

class Node
{
public:
	//force tree
	FLOAT center_of_mass[4];		//center of mass
	FLOAT center_geometric[3];	//geometric center of the node's bounding box
	FLOAT mass;				//total mass in the node
	FLOAT l;					//extension
	
	//tree algorithm
	BBox bounding_box;	//bounding box of the node (min max for each dimension)
	UINT level; 		//level of the node in the tree

	UINT size; 			// size of tree (in nodes), including the subtree underneath it
	UINT address;		// offset in final tree array
	
	std::vector<Particle *> particles;
	std::vector<Chunk> chunks;
	
	Node *left_child;
	Node *right_child;
	
	UINT split_dim;		// dimension along node is split in child nodes. Accelerates sorting into childnodes

	Node() 
	{
		left_child = NULL;
		right_child = NULL;
	};
/*
	Node(const Node& other): mass(other.mass), l(other.l),
			bounding_box(other.bounding_box), level(other.level), size(other.size), address(other.address), particles(other.particles), chunks(other.chunks),
			left_child(other.left_child), right_child(other.right_child), split_dim(split_dim) {
		memcpy(center_of_mass, other.center_of_mass, 3 * sizeof(FLOAT));
		memcpy(center_geometric, other.center_geometric, 3 * sizeof(FLOAT));
	}
*/
};

class Tree
{
friend class TreeWalker;
public:
//private:
	Node *root;
	std::vector<Node> nodelist;
	
	void processLargeNodes(std::vector<Node*> &activelist, std::vector<Node*> &smalllist, std::vector<Node*> &nextlist,
			std::vector<Node> &nodelist, const UINT level);
	void preprocessSmallNodes(std::vector<Node*> &smalllist);
	void processSmallNodes(std::vector<Node*> &activelist, std::vector<Node*> &nextlist, std::vector<Node> &nodelist, UINT& maxLevel);
	Node* preorderTraversal(std::vector<Node> &nodelist, const UINT maxLevel);
	
	void upPass(std::vector<Node*> &activelist);
	void downPass(std::vector<Node*> &activelist, Node* tree);
public:
	Node* buildTree(std::vector<Particle *> &particles, BBox bounding_box);
	
};


#endif
