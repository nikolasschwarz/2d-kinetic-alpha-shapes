#pragma once
#include <cassert>
#include <concepts>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace kinDS
{

template<std::totally_ordered K, typename V>
class AVLTree
{
private:
    enum ChildPosition
    {
        LEFT,
        RIGHT
    };

    class AVLNode
    {
    public:
        // Constructor
        AVLNode(const K& key, const V& value)
            : key(key)
            , value(value)
            , height(1)
        {
        }
        // Destructor
        ~AVLNode()
        {
            if (left)
            {
                delete left;
            }
            if (right)
            {
                delete right;
            }
        }
        // Getters
        K getKey() const { return key; }
        V getValue() const { return value; }
        int getHeight() const { return height; }
        // Setters
        void setValue(const V& value) { this->value = value; }
        void setHeight(int height) { this->height = height; }

        // TODO: rather make iterative
        AVLNode* findNode(const K& key)
        {
            if (key == this->key)
            {
                return this;
            }
            else if (key < this->key && left)
            {
                return left->findNode(key);
            }
            else if (key > this->key && right)
            {
                return right->findNode(key);
            }
            return nullptr;
        }

        AVLNode* findNext()
        {
            // Find the next node in the in-order traversal
            if (right)
            {
                AVLNode* current = right;
                while (current->left)
                {
                    current = current->left;
                }
                return current;
            }
            AVLNode* current = this;
            while (current->parent && current == current->parent->right)
            {
                current = current->parent;
            }
            return current->parent; // NOTE: this can be a nullptr if we are at the rightmost node as current will then be the root node
        }

        // TODO: rather make iterative
        std::pair<AVLNode*, ChildPosition> findInsertionPoint(const K& key)
        {
            if (key < this->key)
            {
                if (!left)
                {
                    return std::make_pair<>(this, ChildPosition::LEFT);
                }
                else
                {
                    return left->findInsertionPoint(key);
                }
            }
            // duplicate keys are inserted to the right
            else
            {
                if (!right)
                {
                    return std::make_pair<>(this, ChildPosition::RIGHT);
                }
                else
                {
                    return right->findInsertionPoint(key);
                }
            }
        }

        void insertChild(AVLNode* child, ChildPosition position)
        {
            assert(child);
            if (position == ChildPosition::LEFT)
            {
                assert(!left); // prevent memory leak
                left = child;
            }
            else
            {
                assert(!right); // prevent memory leak
                right = child;
            }
            child->parent = this;
            child->height = height + 1;
        }

        // call this to remove the node from the tree after all pointers are set
        // delete should not be called as it will also delete the children
        void remove()
        {
            left = nullptr;
            right = nullptr;
            delete this;
        }

        void print()
        {
            if (left)
            {
                left->print();
            }
            std::cout << key << " ";
            if (right)
            {
                right->print();
            }
        }

        inline bool isRoot()
        {
            return parent == nullptr;
        }

        inline bool isLeaf()
        {
            return left == nullptr && right == nullptr;
        }

        K key;
        V value;
        int height;
        AVLNode* left = nullptr;
        AVLNode* right = nullptr;
        AVLNode* parent = nullptr;
    };

    class PrintHelper
    {
    public:
        static int maxDepth(AVLNode* root)
        {
            if (!root)
                return 0;
            return std::max(maxDepth(root->left), maxDepth(root->right)) + 1;
        }

        // Recursive function to fill the 2D vector with node values
        static void fill(AVLNode* root, std::vector<std::vector<std::string>>& res, int row, int col, int depth, int offset)
        {
            if (!root)
                return;
            res[row][col] = std::to_string(root->value);
            int childOffset = offset / 2;

            if (root->left)
            {
                res[row + 1][col - childOffset] = "/";
                fill(root->left, res, row + 2, col - offset, depth, childOffset);
            }
            if (root->right)
            {
                res[row + 1][col + childOffset] = "\\";
                fill(root->right, res, row + 2, col + offset, depth, childOffset);
            }
        }

        // Print function
        static void printTree(AVLNode* root)
        {
            if (!root)
                return;
            int depth = maxDepth(root);
            int height = depth * 2 - 1;
            int width = (1 << depth) * 2 - 1;
            std::vector<std::vector<std::string>> res(height, std::vector<std::string>(width, " "));

            fill(root, res, 0, width / 2, depth, (1 << (depth - 1)));

            for (const auto& row : res)
            {
                for (const auto& val : row)
                {
                    std::cout << val;
                }
                std::cout << "\n";
            }
        }
    };

    AVLNode* root = nullptr;

public:
    // Constructor
    AVLTree() = default;
    // Destructor
    ~AVLTree()
    {
        // Clean up the tree
        if (root)
        {
            delete root;
            root = nullptr;
        }
    }
    // Copy constructor
    AVLTree(const AVLTree& other) = default;
    // Move constructor
    AVLTree(AVLTree&& other) noexcept = default;
    // Copy assignment operator
    AVLTree& operator=(const AVLTree& other) = default;
    // Move assignment operator
    AVLTree& operator=(AVLTree&& other) noexcept = default;

    void insert(const K& key, const V& value)
    {
        // Insert the key-value pair into the tree
        if (!root)
        {
            root = new AVLNode(key, value);
            return;
        }

        std::pair<AVLNode*, ChildPosition> insertionPoint = root->findInsertionPoint(key);

        AVLNode* newNode = new AVLNode(key, value);
        insertionPoint.first->insertChild(newNode, insertionPoint.second);

        // Balance the tree if necessary
    }

    // TODO: Will not work for duplicate keys
    bool exists(const K& key, V& value) const
    {
        // Find the value associated with the key
        // Return true if found, false otherwise
        AVLNode* node = root->findNode(key);
        if (node)
        {
            if (value == node->getValue())
            {
                return true;
            }
        }
        return false;
    }

    std::optional<V> find(const K& key) const
    {
        // Find the value associated with the key
        AVLNode* node = root->findNode(key);
        if (node)
        {
            return node->getValue();
        }

        // Return std::nullopt if not found
        return std::nullopt;
    }

    void remove(const K& key)
    {
        // Remove the key-value pair from the tree
        if (!root)
        {
            return;
        }

        AVLNode* node = root->findNode(key);
        if (!node)
        {
            return; // Key not found
        }

        // If the node is a leaf, simply delete it
        if (node->isLeaf())
        {
            if (node->isRoot())
            {
                root = nullptr;
            }
            else
            {
                if (node->parent->left == node)
                {
                    node->parent->left = nullptr;
                }
                else
                {
                    node->parent->right = nullptr;
                }
            }
        }
        // If the node has one child, replace it with its child
        else if (node->left && !node->right)
        {
            if (node->isRoot())
            {
                root = node->left;
                root->parent = nullptr;
            }
            else
            {
                if (node->parent->left == node)
                {
                    node->parent->left = node->left;
                }
                else
                {
                    node->parent->right = node->left;
                }
                node->left->parent = node->parent;
            }
        }
        else if (!node->left && node->right)
        {
            if (node->isRoot())
            {
                root = node->right;
                root->parent = nullptr;
            }
            else
            {
                if (node->parent->left == node)
                {
                    node->parent->left = node->right;
                }
                else
                {
                    node->parent->right = node->right;
                }
                node->right->parent = node->parent;
            }
        }
        // If the node has two children, find the in-order successor
        else
        {
            AVLNode* Y = node->findNext();
            if (Y->parent != node)
            {
                // replace Y by its right child
                AVLNode* X = Y->right;
                if (X)
                {
                    X->parent = Y->parent;
                }

                // In this case also reassign the right child
                Y->right = node->right;
            }
            // Y cannot have a left child, replace node with Y and inherit the left child from node
            if (node->isRoot())
            {
                root = Y;
            }

            Y->parent = node->parent;
            Y->left = node->left;
        }

        // free memory from the node without deleting the (now invalid) children
        node->remove();
        // Update the height of the node
        // TODO
        // Balance the tree if necessary
        // TODO
    }

    void print() const
    {
        // Print the tree in-order
        if (root)
        {
            root->print();
        }
    }

    void printTreeStructure() const
    {
        // Print the tree structure
        if (root)
        {
            PrintHelper::printTree(root);
        }
    }

private:
    int height(AVLNode* node) const
    {
        return node ? node->height : 0;
    }
    int balanceFactor(AVLNode* node) const
    {
        return height(node->left) - height(node->right);
    }
    void rotateLeft(AVLNode*& node)
    {
        // Perform left rotation
    }
    void rotateRight(AVLNode*& node)
    {
        // Perform right rotation
    }
    void balance(AVLNode*& node)
    {
        // Balance the tree
    }
};
} // namespace kinDS
