#include<iostream>
#include"displayTets.h"

int main(int argc, char * argv[]){
    std::string Node(argv[2]), Face(argv[1]);
    displayTet(Node.c_str(),Face.c_str());
    return 0;
}
