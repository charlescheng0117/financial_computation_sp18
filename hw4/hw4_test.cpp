

void testMaxStList() {
	vector<double> st = {1, 2, 3, 4, 5, -5.5};

	MaxStList s = MaxStList(st);

	s.display();

	cout << "max st: " << s.max_St << "\n";
	cout << "size: " << s.size << "\n";

	s.add_St(100);

	s.display();

	cout << "max st: " << s.max_St << "\n";
	cout << "size: " << s.size << "\n";

	s.add_St(2);

	s.display();

	cout << "max st: " << s.max_St << "\n";
	cout << "size: " << s.size << "\n";

}

void testNode() {
	long double r = 0.1, sigma = 0.25, q = 0, u = 1.1224, d = 0.8909, p = 0.5073;
	double t = 0, T = 0.25, dT = 0.083333;
	int n = 3;

	double St = 50;

	cout << "root node.";
	Node root = Node(St);
	root.display();

	cout << "test node:";
	Node test = Node();
	test.display();

	cout << "size should be 0. test is: " << test.S_max_vals->size;

	Node test_cons = *(new Node(100));

	test_cons.display();

	return;
}

void testInherit() {
	vector<double> maxs_1 = {1, 4, 7, 8, 9, 30};
	Node n1 = Node(1, maxs_1);

	n1.display();

	vector<double> maxs_2 = {1, 3, 4, 10, 25.5, 30, 60};

	Node n2;

	n2.display();

	n2.inheritSmax(&n1);
	n2.display();

	Node n3 = Node(25, maxs_2);

	cout << "adklsjdfakl; NODE3\n";
	n3.display();

	n3.inheritSmax(&n1);

	n3.display();

}

void testSpawnChild() {
	long double r = 0.1, sigma = 0.25, q = 0, u = 1.1224, d = 0.8909, p = 0.5073;
	double t = 0, T = 0.25, dT = 0.083333;
	int n = 3;

	double St = 50;

	printf("u: %Lf, d: %Lf, p: %Lf, dT: %f\n", u, d, p, dT);

	vector<double>* root_S_max = new vector<double>;
	root_S_max->push_back(St);

	Node root = Node(St, *root_S_max);

	cout << "\nLevel 0\n";
	root.display();

	// constructing level 1
	root.addUpChild(u);
	root.addDownChild(d);

	Node* child_1_0 = root.down_child;
	Node* child_1_1 = root.up_child;
	child_1_0->up_sib = child_1_1;

	cout << "\nLevel 1\n";
	child_1_1->display();
	child_1_0->display();

	// constructing level 2
	child_1_0->addDownChild(d);
	child_1_0->addUpChild(u);     // also the down child of child_1_1

	child_1_0->up_child->inheritSmax(child_1_1);
	child_1_1->down_child = child_1_0->up_child;

	child_1_1->addUpChild(u);
	
	Node* child_2_2 = child_1_1->up_child;
	Node* child_2_1 = child_1_0->up_child;
	Node* child_2_0 = child_1_0->down_child;

	child_2_0->up_sib = child_2_1;
	child_2_1->up_sib = child_2_2;

	cout << "\nLevel 2\n";
	child_2_2->display();
	child_2_1->display();
	child_2_0->display();

	// constructing level 3

	child_2_0->addDownChild(d);
	child_2_0->addUpChild(u);

	child_2_0->up_child->inheritSmax(child_2_1);
	child_2_1->down_child = child_2_0->up_child;

	child_2_1->addUpChild(u);
	child_2_1->up_child->inheritSmax(child_2_2);
	child_2_2->down_child = child_2_1->up_child;

	child_2_2->addUpChild(u);

	Node* child_3_3 = child_2_2->up_child;
	Node* child_3_2 = child_2_1->up_child;
	Node* child_3_1 = child_2_0->up_child;
	Node* child_3_0 = child_2_0->down_child;

	child_3_0->up_sib = child_3_1;
	child_3_1->up_sib = child_3_2;
	child_3_2->up_sib = child_3_3;

	cout << "\nLevel 3\n";
	child_3_3->display();
	child_3_2->display();
	child_3_1->display();
	child_3_0->display();

	Node* bottom = goBottomDown(&root);
	Node* cur = bottom;
	Node* next_bottom = cur->up_parent;

	while (cur != NULL) {
		cur->calcOptionVal(u, r, q, dT, p);
		cur->display();

		cur = cur->up_sib;
	}

	cur = next_bottom;
	next_bottom = cur->up_parent;

	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	cur = cur->up_sib;
	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	cur = cur->up_sib;
	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	cur = next_bottom;
	next_bottom = cur->up_parent;

	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	cur = cur->up_sib;
	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	cur = next_bottom;
	cur->calcOptionVal(u, r, q, dT, p);
	cur->display();

	deleteTree(&root);
}

void testRound() {
	cout << roundDouble(0.122401, 4) << '\n';
	cout << roundDouble(0.890947, 4) << '\n';
	cout << roundDouble(0.507319, 4) << '\n';
	cout << roundDouble(0.083333, 4) << '\n';
	cout << roundDouble(44.547363, 2) << '\n';
	cout << roundDouble(39.689350, 2) << '\n';
	cout << roundDouble(56.120045, 2) << '\n';
	cout << roundDouble(62.989189, 2) << '\n';
}
