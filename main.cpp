#include <iostream>
#include "kinDS/AVLTree.hpp"

int main()
{
    kinDS::AVLTree<int, int> test;

    test.insert(1, 1);
    test.insert(3, 3);
    test.insert(4, 4);
    test.insert(5, 5);
    test.insert(2, 2);
    
    test.printTreeStructure();
}
