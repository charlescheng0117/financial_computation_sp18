#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <queue>

#define NODE_CREATE 0
#define NODE_DESTROYED 0

using namespace std;

int NODE_TOTAL = 0;

struct pair_option {
	double euro;
	double amer;

    pair_option() {}

    pair_option(double e, double a) {
        euro = e;
        amer = a;
    }
};

double roundDouble(double val, int d) {
	return round( val * pow(10.0, d) ) / pow(10.0, d);
}

void displayVec(vector<double>& v) {
	for (int i = 0; i < v.size(); i += 1) {
		cout << v[i] << " ";
	}
	cout << "\n";
}

class MaxStList {
public:
	vector<double> St_list;
	double max_St;
	int list_size;

	MaxStList() {
		max_St = -1;
		list_size = 0;
	}

	MaxStList (vector<double>& lst) {
		list_size = lst.size();
        St_list = vector<double>(list_size);

		for (int i = 0; i < list_size; i += 1) {
		    St_list[i] = lst[i];
        }
		//cout << "size is: " << size << '\n';
		max_St = *max_element(St_list.begin(), St_list.end());
	}

	void add_St(double St) {
		if (St > max_St) {
			St_list.push_back(St);
			max_St = St;
			list_size += 1;
		}
	}

	void display() {
		if (list_size == 0) {
			cout << "empty\n";
			return;
		}

		for (int i = 0; i < list_size; i += 1) {
			cout << St_list[i] << " ";
		}
		cout << '\n';
	}
};


class Node {
public:
	double St = -1;             		  // Current price
	int depth = 0;              		  // Depth of this Node in the tree
	bool isLeaf = true;

	MaxStList S_max_vals;        		  // Record each possible S_max
	vector<double> option_european;  	  // Value of the european option
	vector<double> option_american;      // Value of the american option

	Node *up_parent;
	Node *down_parent;
	Node *up_sib;
	Node *down_sib;
	Node *up_child;
	Node *down_child;

	Node () {
		St = -1, depth = 0, isLeaf = true;

		up_parent = down_parent = NULL;
		up_sib = down_sib = NULL;
		up_child = down_child = NULL;

		NODE_TOTAL += 1;
	}

	Node (double St) { // root constructor, maybe
		this->St = St, depth = 0, isLeaf = false;

		NODE_TOTAL += 1;
	}

    void create_Node(double new_St) {
        St = new_St;
        depth = 0;
        isLeaf = false;
    }

    void create_Node(double new_St, vector<double>& St_list) {
		St = new_St;	
		depth = 0;
		isLeaf = false;

		// Recording possible Smax
		int size_list = St_list.size();
		bool isStAdded = false;

		for (int i = 0; i < size_list; i += 1) {
			double S_max = St_list[i];
			if (S_max >= St) {
				S_max_vals.add_St(S_max);
			} else { // parent_S_max < St
				if (!isStAdded)
					S_max_vals.add_St(St);
					isStAdded = true;
			}
		}

		up_parent = down_parent = NULL;
		up_sib = down_sib = NULL;
		up_child = down_child = NULL;

		NODE_TOTAL += 1;
    
    }

	Node (double St, vector<double>& St_list) { // constructor for root node
		this->St = St;	
		depth = 0;
		isLeaf = false;

		// Recording possible Smax
		int size_list = St_list.size();
		bool isStAdded = false;

		for (int i = 0; i < size_list; i += 1) {
			double S_max = St_list[i];
			if (S_max >= St) {
				//cout << "S_max now is: " << S_max << "\n";
				S_max_vals.add_St(S_max);
			} else { // parent_S_max < St
				if (!isStAdded)
					S_max_vals.add_St(St);
					isStAdded = true;
			}
		}

		up_parent = down_parent = NULL;
		up_sib = down_sib = NULL;
		up_child = down_child = NULL;

		NODE_TOTAL += 1;
	}

	~Node() { 
		//delete option_european;
		//delete option_american;
		//cout << "A Node is destroyed.\n";
		NODE_TOTAL -= 1;
	}


	void inheritSmax(Node* parent) {
		if (parent == NULL) { 					// the up-most node
			cout << "No parent node here.\n";
		}
		if (S_max_vals.list_size == 0) { 			// Just consider the case listed in Alg (i)
			vector<double> parent_S_max_list = (parent->S_max_vals.St_list);
			int num_parent = parent_S_max_list.size();
			bool isStAdded = false;

			for (int i = 0; i < num_parent; i += 1) {
				double parent_S_max = parent_S_max_list[i];

				if (parent_S_max >= St) {
					S_max_vals.add_St(parent_S_max);
				} else { // parent_S_max < St
					if (!isStAdded) {
						S_max_vals.add_St(St);
						isStAdded = true;
					}
				}
			}
		} else { // have to combine two Smax list.
			vector<double> my_S_max_list = (S_max_vals.St_list);
			vector<double> parent_S_max_list = (parent->S_max_vals.St_list);
			vector<double> combined_S_max_list;

			sort(my_S_max_list.begin(), my_S_max_list.end());
			sort(parent_S_max_list.begin(), parent_S_max_list.end());

			int my_size = my_S_max_list.size(), parent_size = parent_S_max_list.size();
			int i = 0, j = 0;

			while (i < my_size && j < parent_size) {
				if (my_S_max_list[i] > parent_S_max_list[j]) {

					if (parent_S_max_list[j] >= this->St) // only inherit from parents if Smax >= St
						combined_S_max_list.push_back(parent_S_max_list[j]);
					j += 1;

				} else if (my_S_max_list[i] < parent_S_max_list[j]) {

					if (my_S_max_list[i] >= this->St) // only inherit from parents if Smax >= St
						combined_S_max_list.push_back(my_S_max_list[i]);
					i += 1;

				} else { // equal

					if (parent_S_max_list[j] >= this->St) // only inherit from parents if Smax >= St
						combined_S_max_list.push_back(parent_S_max_list[j]);
					i += 1;
					j += 1;
				}
			}

			if (i == my_size) {
				while (j < parent_size) {
					if (parent_S_max_list[j] >= this->St) // only inherit from parents if Smax >= St
						combined_S_max_list.push_back(parent_S_max_list[j]);
					j += 1;
				}
			} else {
				while (i < my_size) {
					if (my_S_max_list[i] >= this->St) // only inherit from parents if Smax >= St
						combined_S_max_list.push_back(my_S_max_list[i]);
					i += 1;
				}
			}
			// update size, max_list
			S_max_vals = MaxStList(combined_S_max_list);

			return;
		}
	}

	void addUpChild(double u) {
		//double uS_t = roundDouble(u * St, 2);
		double uS_t = u * St;
        //printf("uS_t = %f\n", uS_t);
		//up_child = new Node(uS_t);       // create a upper child.
        up_child = new Node;
        up_child->create_Node(uS_t);
		up_child->depth = depth + 1;     // update depth
		up_child->down_parent = this;    // for this child, `this' is its `down_parent'
		up_child->inheritSmax(this);     // inherit Smax from `this'
	}

	void addDownChild(double d) {
		//double dS_t = roundDouble(d * St, 2);
		double dS_t = d * St;
		//down_child = new Node(dS_t);     // create a lower child.
        down_child = new Node;
        down_child->create_Node(dS_t);
		down_child->depth = depth + 1;   // update depth
		down_child->up_parent = this;    // for this child, `this' is its `up_parent'
		down_child->inheritSmax(this);   // inherit Smax from `this'
	}

	pair_option findOptionVal(double St, double S_max, double u) { 
	// given St, find the corresponding index in S_max_list. Then find the option val.
		vector<double> vec_S_max = (S_max_vals.St_list);
		int idx = 0;
		bool found = false;
		//printf("St = %f\nvec_S_max: \n", St);
		//displayVec(vec_S_max);

		for (int i = 0; i < vec_S_max.size(); i += 1) {
			//if (St == vec_S_max[i]) {
			//if (abs(S_max - vec_S_max[i]) <= 0.01) {
			if (abs(S_max - vec_S_max[i]) <= 0.00001) {
				idx = i;
				found = true;
				break;
			}
		}

		pair_option ret;
		if(found) {
			ret.euro = option_european[idx];
			ret.amer = option_american[idx];
			return ret;
		} else {   // find the put value for S_max = S_t * u
			
			//double uSt = roundDouble(St * u, 2);
			double uSt = St * u;

			//return this->findOptionVal(St, uSt, u);
            return findOptionVal(St, uSt, u);
		}
	}

	void calcOptionVal(double u, double r, double q, double dt, double p) {
		// if leaf:
		vector<double> vec_S_max = (S_max_vals.St_list);
		if (isLeaf) {
			for (int i = 0; i < vec_S_max.size(); i += 1) { // With S_max_i, calculate the put value for put_i
				double S_max = vec_S_max[i];
				option_european.push_back(max(S_max - St, 0.0));
				option_american.push_back(max(S_max - St, 0.0));
			}
		} else {
			Node* _up_child = this->up_child;
			Node* _down_child = this->down_child;

			for (int i = 0; i < vec_S_max.size(); i += 1) { // With S_max_i, calculate the put value for put_i
				double S_max = vec_S_max[i];

				//pair_option* upOptVal = _up_child->findOptionVal(St, S_max, u);     // find up-pair of option vals
				pair_option upOptVal = _up_child->findOptionVal(St, S_max, u);     // find up-pair of option vals
				//pair_option* downOptVal = _down_child->findOptionVal(S_max, u); // find down-pair of opiton vals
				pair_option downOptVal = _down_child->findOptionVal(St, S_max, u); // find down-pair of opiton vals

				// calculate european option value
				double eu_opt_val = exp(-r * dt) * (p * (upOptVal.euro) + (1 - p) * (downOptVal.euro) );
				//eu_opt_val = roundDouble(eu_opt_val, 2); 
				option_european.push_back(eu_opt_val);

				// calculate american option value
				double am_opt_val = exp(-r * dt) * (p * (upOptVal.amer) + (1 - p) * (downOptVal.amer) );
				double EV = S_max - St;
				//am_opt_val = roundDouble(am_opt_val, 2);

				if (EV > am_opt_val) {				       // early excercise
					option_american.push_back(EV);
				} else {
					option_american.push_back(am_opt_val);
				}

			}
		}
	}

	void display() {
		printf("\nNode info:\n");
		printf("St: %f Depth: %d IsLeaf: %d\n", St, depth, isLeaf);
	    printf("S_max values are: ");
		if (S_max_vals.list_size != 0) {
			S_max_vals.display();
			if (option_european.size() != 0) {
				printf("European Option vals are:  ");
				displayVec(option_european);
				printf("American Option vals are:  ");
				displayVec(option_american);
			}
		} else {
			cout << "S_max's are empty. \n";
		}
		//displayVec(*option_european);
	}
};

void displayTree(Node* root, int n) {   // n: total depth
	queue<Node*> nodeQueue;
	nodeQueue.push(root);
	int cur_depth = root->depth;
	cout << "\nDisplaying Tree.\n";
	while (cur_depth <= n) {                       // in depth n, there are (n + 1) nodes to spawn a child.
		printf("\nlevel %d:\n", cur_depth);
		Node* bottom_node = nodeQueue.front();    // the lower most node to spawn 2 child.
		nodeQueue.pop();
		bottom_node->display();

		nodeQueue.push(bottom_node->down_child);
		nodeQueue.push(bottom_node->up_child);

		Node* prev_node = bottom_node->up_child;  // `prev_node' should inherit from the next node
		 										  // popped from the queue.

		for (int i = 1; i <= cur_depth; i += 1) {
			Node* next_node = nodeQueue.front();  // pop out the next node to spawn up child.
			nodeQueue.pop();
			next_node->display();

			nodeQueue.push(next_node->up_child);
			prev_node = next_node->up_child;      // this node becomes the `prev_node' 
		}
		cur_depth += 1;
	}
}

Node* goBottomDown(Node* root) {  // Return a pointer to the ``bottom-most'' node in the tree.
	Node* cur_Node = root;
	Node* bottom_most;
	//Node* prev_Node;

	while (cur_Node->down_child != NULL) { // Walk down through the pointer
		cur_Node = cur_Node->down_child;
	}
	bottom_most = cur_Node;

	while (cur_Node->up_sib != NULL) { // Set up the isLeaf variable. 
		cur_Node->isLeaf = true;
		cur_Node = cur_Node->up_sib;
	}
	cur_Node->isLeaf = true;

	return bottom_most;
}

void deleteTreeHelper(Node* head) {    // Delete the entire links among the siblings.
	if (head->up_sib == NULL) {
		delete head;
		return;
	}
	deleteTreeHelper(head->up_sib);
	delete head;
	return;
}

void deleteTree(Node* root) {
	Node* bottom = goBottomDown(root);
	Node* cur = bottom;
	Node* next_bottom;

	while (cur != root) {
		next_bottom = cur->up_parent;
		deleteTreeHelper(cur);
		cur = next_bottom;
	}
	return;
}

pair_option binomLookback( double St, double r, double q, double sigma, double t, double T,
				  double S_max_t, int n, int simulations, int repetitions ) {
	printf("Number of periods: %d\n", n);
	printf("St = %f, r = %f, q = %f, sigma = %f, t = %f, T = %f, S_max_t = %f\n", St, r, q, sigma, t, T, S_max_t);

	double u, d, p, dT;

	dT = (T - t) / n;
	u = exp(sigma * sqrt(dT));
	d = 1/u;
	p = (exp((r-q) * dT) - d)/(u - d);

	// round up
	//dT = roundDouble(dT, 4), u = roundDouble(u, 4), d = roundDouble(d, 4), p = roundDouble(p, 4);
	printf("u: %f, d: %f, p: %f, dT: %f\n", u, d, p, dT);

	vector<double> root_S_max;
	root_S_max.push_back(S_max_t);

	Node root = Node(St, root_S_max);

	//root.display();

	queue<Node*> nodeQueue; // a FIFO queue
	nodeQueue.push(&root);
	int cur_depth = root.depth;
	
	// TODO: SET UP THE CORRECT LINK BETWEEN YOUR CHILDREN.

	while (cur_depth < n) {                       				 // in depth n, there are (n + 1) nodes to spawn a child.
		Node* bottom_node = nodeQueue.front();    				 // the lower most node to spawn 2 child.
		nodeQueue.pop();

		bottom_node->addDownChild(d);			  				 // spawn the down-child
		nodeQueue.push(bottom_node->down_child);

		bottom_node->addUpChild(u);				  				 // spawn the up-child
		nodeQueue.push(bottom_node->up_child);

		bottom_node->down_child->up_sib = bottom_node->up_child; // setup the up-sibling link.

		Node* prev_up_child = bottom_node->up_child;  		     // `prev_node'(up_child of `this') should inherit 
																 // from the next node(up_sib of `this')
		 										  				 // popped from the queue.
		for (int i = 1; i <= cur_depth; i += 1) {
			Node* next_up_sib = nodeQueue.front();  			 // pop out the next node to spawn up child.
			nodeQueue.pop();

			prev_up_child->inheritSmax(next_up_sib);    		 // inheriting
			next_up_sib->down_child = prev_up_child;             // setup the child link
			prev_up_child->up_parent = next_up_sib;				 // setup the parent link

			next_up_sib->addUpChild(u);             				 // spawn the up-child
			prev_up_child->up_sib = next_up_sib->up_child;		     // setup the up-sibling link.

			prev_up_child = next_up_sib->up_child;      			 // this node becomes the `prev_node'  
			nodeQueue.push(prev_up_child);             				 // should push into the queue.
		}
		cur_depth += 1;
	}
	Node* bottom = goBottomDown(&root);

	// the sub-routine to calculate the option price.
	Node* cur_Node = bottom;
	Node* next_bottom;

	// for delete node
	Node* top_Node;
	Node* cur_bottom;
	
	while (cur_Node->depth != 0) {
		cur_bottom  = cur_Node;
		next_bottom = cur_Node->up_parent;

		int i = 0;
		while (cur_Node != NULL) {
			cur_Node->calcOptionVal(u, r, q, dT, p);
            
            // debug
            //printf("#depth %d: ", i);
            //cur_Node->display();

			/*
			printf("i = %d\n", i);
			if (cur_Node->up_sib == NULL && cur_Node->isLeaf == false) {
				printf("Start deleting the link");
				top_Node = cur_Node;
			}
			printf("i = %d\n", i);
			*/
			cur_Node = cur_Node->up_sib;
			i += 1;
		}
		//displayTree(&root, n);

		/*
		printf("Start deleting the link");
		if (cur_bottom != NULL)
			deleteTreeHelper(cur_bottom->down_child);
		*/

		cur_Node = next_bottom;
	}

	// calculate the option value of root
	cur_Node->calcOptionVal(u, r, q, dT, p);

	//displayTree(&root, n);

	pair_option option_val;
	option_val.euro = (root.option_european)[0]; // european option
	option_val.amer = (root.option_american)[0]; // american option

	//cout << "Euro: " << option_val.euro << "\n";
	//cout << "Amer: " << option_val.amer << "\n";
	//deleteTree(&root);
	return option_val;
}

int main(int argc, char const **argv) {
	// load data
	if (argc != 2) {
        printf("Usage: ./binom <input_file>\n");
		return -1;
	}

	const char* in_file = argv[1];
	ifstream in(in_file);

	double St, r, q, sigma, t, T, S_max_t;
	int n, simulations, repetitions;

	
	in >> St >> r >> q >> sigma >> t >> T >> S_max_t;
	in >> n >> simulations >> repetitions;


	/*
	St = 50, r = 0.1, q = 0.05, sigma = 0.4, t = 1, T = 4, S_max_t = 60;
	n = 100, simulations = 1000, repetitions = 20;*/
    
	printf("Basic requirement: Binomial Tree Model\n");
	pair_option option_val = binomLookback(St, r, q, sigma, t, T, S_max_t, n, simulations, repetitions);
	printf("### Answer for Binomial Tree Model ###\n");
    printf("European: %f\nAmerican: %f\n", option_val.euro, option_val.amer);

	return 0;
}
