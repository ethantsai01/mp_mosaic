/**
 * @file kdtree.cpp
 * Implementation of KDTree class.
 */

#include <utility>
#include <algorithm>

using namespace std;

template <int Dim>
bool KDTree<Dim>::smallerDimVal(const Point<Dim>& first,
                                const Point<Dim>& second, int curDim) const
{
    /**
     * @todo Implement this function!
     */
     if(first[curDim] < second[curDim]) {
       return true;
     }
     else if(first[curDim] > second[curDim]) {
       return false;
     }
    return first < second;
}

template <int Dim>
bool KDTree<Dim>::shouldReplace(const Point<Dim>& target,
                                const Point<Dim>& currentBest,
                                const Point<Dim>& potential) const
{
    /**
     * @todo Implement this function!
     */
     double currDist = distance(target, currentBest); //distance of currentBest to target
     double checkDist = distance(target, potential); //distance of potential to target

     if(checkDist < currDist) {
       return true;
     }
     else if(checkDist > currDist) {
       return false;
     }
     return checkDist < currDist;
}

template <int Dim>
int KDTree<Dim>::distance(const Point<Dim>& target, const Point<Dim>& input) const {
    int distance = 0;
    for(int i = 0; i < Dim; i++) {
      distance += (target[i] - input[i])*(target[i] - input[i]);
    }
    return distance;
}

template <int Dim>
KDTree<Dim>::KDTree(const vector<Point<Dim>>& newPoints)
{
    /**
     * @todo Implement this function!
     */
     if(newPoints.empty() == true) {
       root = NULL;
     }
     else{
       vector<Point<Dim>> tempPoints;
       for(unsigned long i = 0; i < newPoints.size(); i++) {
         tempPoints.push_back(newPoints[i]);
       }
       root = build(tempPoints, 0, tempPoints.size() - 1, 0);
     }
}

template<int Dim>
typename KDTree<Dim>::KDTreeNode* KDTree<Dim>::build(vector<Point<Dim>>& newPoints, int left, int right, int dimension) {
  if(left > right) {
    return NULL;
  }

  int med = (left + right)/2;
  Point<Dim> median = pick(newPoints, left, right, med, dimension);
  KDTreeNode* subRoot = new KDTreeNode(median);

  size += 1;
  subRoot -> left = build(newPoints, left, med - 1, (dimension + 1)%Dim);
  subRoot -> right = build(newPoints, med + 1, right, (dimension + 1)%Dim);
  return subRoot;
}

template<int Dim>
Point<Dim> KDTree<Dim>::pick(vector<Point<Dim>>& newPoints, int left, int right, int k, int dimension) {
  if(left == right) {
    return newPoints[left];
  }

  int pivot = (left + right)/2;
  pivot = part(newPoints, left, right, pivot, dimension);
  if(k == pivot) {
    return newPoints[k];
  }
  else if(k < pivot) {
    return pick(newPoints, left, pivot - 1, k, dimension);
  }
  else {
    return pick(newPoints, pivot + 1, right, k, dimension);
  }
}

template<int Dim>
int KDTree<Dim>::part(vector<Point<Dim>>& newPoints, int left, int right, int pivot, int dimension) {
  Point<Dim> pivotVal = newPoints[pivot];
  swap(newPoints, pivot, right);
  int store = left;
  for(int i = left; i < right; i++) {
    if(smallerDimVal(newPoints[i], pivotVal, dimension)) {
      swap(newPoints, store, i);
      store++;
    }
  }
  swap(newPoints, right, store);
  return store;
}

template <int Dim>
void KDTree<Dim>::swap(vector<Point<Dim>>& newPoints, int left, int right) {
  Point<Dim> temp = newPoints[left];
  newPoints[left] = newPoints[right];
  newPoints[right] = temp;
}

template <int Dim>
KDTree<Dim>::KDTree(const KDTree<Dim>& other) {
  /**
   * @todo Implement this function!
   */
   copy(root, other -> root);
   size = other.size;
}

template <int Dim>
void KDTree<Dim>::copy(KDTreeNode*& current, KDTreeNode*& other) {
  if(other == NULL) {
     return;
   }

   current = new KDTreeNode(other -> print);
   copy(current -> left, other -> left);
   copy(current -> right, other -> right);
}

template <int Dim>
const KDTree<Dim>& KDTree<Dim>::operator=(const KDTree<Dim>& rhs) {
  /**
   * @todo Implement this function!
   */
  delete(root);
  copy(root, rhs -> root);
  size = rhs.size;
  return *this;
}

template <int Dim>
KDTree<Dim>::~KDTree() {
  /**
   * @todo Implement this function!
   */
   delete(root);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighbor(const Point<Dim>& query) const
{
    /**
     * @todo Implement this function!
     */

    return findNearestNeighborHelper(query, root, 0);
}

template <int Dim>
Point<Dim> KDTree<Dim>::findNearestNeighborHelper(const Point<Dim>& query, typename KDTree<Dim>::KDTreeNode* subRoot, int dimension) const {
  if(subRoot -> left == NULL & subRoot -> right == NULL) {
    return subRoot -> point;
  }

  int nextDim = ((dimension + 1)%Dim);
  Point<Dim> currentBest = subRoot -> point;
  Point<Dim> potBest = currentBest;
  bool direction = smallerDimVal(query, currentBest, dimension);
  if(direction && subRoot -> left == NULL) {
    potBest = findNearestNeighborHelper(query, subRoot -> right, nextDim);
  }
  else if(direction && subRoot -> left != NULL) {
    potBest = findNearestNeighborHelper(query, subRoot -> left, nextDim);
  }
  else if(!direction && subRoot -> right == NULL) {
    potBest = findNearestNeighborHelper(query, subRoot -> left, nextDim);
  }
  else if(!direction && subRoot -> right != NULL) {
    potBest = findNearestNeighborHelper(query, subRoot -> right, nextDim);
  }

  if(shouldReplace(query, currentBest, potBest)) {
    currentBest = potBest;
  }

  int bestDist = distance(currentBest, query); //flipped currentBest and query
  int dimDist = (subRoot -> point[dimension] - query[dimension])*(subRoot -> point[dimension] - query[dimension]);
  if(dimDist <= bestDist) {
    if(!direction && subRoot -> left != NULL) {
      potBest = findNearestNeighborHelper(query, subRoot -> left, nextDim);
      if(shouldReplace(query, currentBest, potBest)) {
        currentBest = potBest;
      }
    }
    else if(direction && subRoot -> right != NULL) {
      potBest = findNearestNeighborHelper(query, subRoot -> right, nextDim);
      if(shouldReplace(query, currentBest, potBest)) {
        currentBest = potBest;
      }
    }
  }
  return currentBest;
}





// bool state;
// Point<Dim> currentBest = subRoot -> point;
// if(subRoot -> left == NULL & subRoot -> right == NULL) {
//   return subRoot -> point;
// }
//
// if(smallerDimVal(query, subRoot -> point, dimension)) {
//   if(subRoot -> left == NULL) {
//     currentBest = findNearestNeighborHelper(query, subRoot -> right, (dimension + 1) % Dim);
//   }
//   if(subRoot -> left != NULL) {
//     currentBest = findNearestNeighborHelper(query, subRoot -> left, (dimension + 1) % Dim);
//   }
//   state = true;
// }
// else {
//   if(subRoot -> right == NULL) {
//     currentBest = findNearestNeighborHelper(query, subRoot -> left, (dimension + 1) % Dim);
//   }
//   if(subRoot -> right != NULL) {
//     currentBest = findNearestNeighborHelper(query, subRoot -> right, (dimension + 1) % Dim);
//   }
//   state = false;
// }
//
// if(shouldReplace(query, currentBest, subRoot -> point)) {
//   currentBest = subRoot -> point;
// }
//
// double radius = 0;
//
// for(int i = 0; i < Dim; i++) {
//   radius += (currentBest[i] = query[i]) * (currentBest[i] = query[i]);
// }
//
// double distance = (subRoot -> point[dimension] - query[dimension]) * (subRoot -> point[dimension] - query[dimension]);
//
// if(distance <= radius) {
//   if(state == true && subRoot -> right != NULL) {
//     Point<Dim> pointCheck = findNearestNeighborHelper(query, subRoot -> right, (dimension + 1) % Dim);
//     if(shouldReplace(query, currentBest, pointCheck)) {
//       currentBest = pointCheck;
//     }
//   }
//   else if(state == false && subRoot -> left != NULL) {
//     Point<Dim> pointCheck = findNearestNeighborHelper(query, subRoot -> left, (dimension + 1) % Dim);
//     if(shouldReplace(query, currentBest, pointCheck)) {
//       currentBest = pointCheck;
//     }
//   }
// }

// return currentBest;
