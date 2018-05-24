#include <float.h>
#include <assert.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"
#include "error_dialog.h"

namespace CS248 {

VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
  // TODO: (meshEdit)
  // This method should split the given edge and return an iterator to the
  // newly inserted vertex. The halfedge of this vertex should point along
  // the edge that was split, rather than the new edges.

  HalfedgeIter h0 = e0->halfedge();
  bool is_boundary0 = h0->isBoundary(), is_boundary1 = h0->twin()->isBoundary();

  // Check triangles
  HalfedgeIter h;
  int count0 = 0;
  if (!is_boundary0) {
    h = h0;
    do {
      count0++;
      h = h->next();
    } while (h != h0);
  }

  int count1 = 0;
  if (!is_boundary1) {
    h = h0->twin();
    do {
      count1++;
      h = h->next();
    } while (h != h0->twin());
  }

  if (is_boundary0 && is_boundary1) {
    showError("Both half-edges are boundary. Something must go wrong.");
    return h0->vertex();
  }

  if ((!is_boundary0 && count0 != 3) || (!is_boundary1 && count1 != 3)) {
    showError("splitEdge() cannot be applied to non-triangle meshes.");
    return h0->vertex();
  }

  HalfedgeIter h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12, h13, h14, h15;
  EdgeIter e1, e2, e3, e4, e5, e6, e7;
  VertexIter v0, v1, v2, v3, v4;
  FaceIter f0, f1, f2, f3;

  // Collect original elements
  if (!is_boundary0) {
    h0 = e0->halfedge();
    h1 = h0->next();
    h2 = h1->next();
    h3 = h0->twin();
    h6 = h1->twin();
    h7 = h2->twin();
    e1 = h1->edge();
    e2 = h2->edge();
    v0 = h0->vertex();
    v1 = h1->vertex();
    v2 = h2->vertex();
  }
  
  if (!is_boundary1) {
    h0 = e0->halfedge();
    h3 = h0->twin();
    h4 = h3->next();
    h5 = h4->next();
    h8 = h4->twin();
    h9 = h5->twin();
    e3 = h4->edge();
    e4 = h5->edge();
    v1 = h3->vertex();
    v0 = h4->vertex();
    v3 = h5->vertex();  
  }

  f0 = h0->face();
  f1 = h3->face();
  
  // Allocate new elements
  v4 = newVertex();
  v4->position = (v0->position + v1->position) / 2.0;
  e5 = newEdge();
  h10 = newHalfedge();
  h11 = newHalfedge();
  if (!is_boundary0) {
    e6 = newEdge();
    h12 = newHalfedge();
    h13 = newHalfedge();
    f3 = newFace();
  }
  if (!is_boundary1) {
    e7 = newEdge();
    h14 = newHalfedge();
    h15 = newHalfedge();
    f2 = newFace();
  }

  // Reassign elements
  if (!is_boundary0) {
    h0->setNeighbors(h1, h3, v4, e0, f0);
    h1->setNeighbors(h12, h6, v1, e1, f0);
    h2->setNeighbors(h10, h7, v2, e2, f3);
    h6->setNeighbors(h6->next(), h1, v2, e1, h6->face());
    h7->setNeighbors(h7->next(), h2, v0, e2, h7->face());
    h10->setNeighbors(h13, h11, v0, e5, f3);
    h12->setNeighbors(h0, h13, v2, e6, f0);
    h13->setNeighbors(h2, h12, v4, e6, f3);
    v2->halfedge() = h2;
    e6->halfedge() = h12;
    f3->halfedge() = h10;
  }

  if (!is_boundary1) {
    h3->setNeighbors(h14, h0, v1, e0, f1);
    h4->setNeighbors(h15, h8, v0, e3, f2);
    h5->setNeighbors(h3, h9, v3, e4, f1);
    h8->setNeighbors(h8->next(), h4, v3, e3, h8->face());
    h9->setNeighbors(h9->next(), h5, v1, e4, h9->face());
    h11->setNeighbors(h4, h10, v4, e5, f2);
    h14->setNeighbors(h5, h15, v4, e7, f1);
    h15->setNeighbors(h11, h14, v3, e7, f2);
    v3->halfedge() = h5;
    e7->halfedge() = h15;
    f2->halfedge() = h11;
  }

  v0->halfedge() = h10;
  v1->halfedge() = h3;
  v4->halfedge() = h0;
  e0->halfedge() = h0;
  e5->halfedge() = h10;
  f0->halfedge() = h0;
  f1->halfedge() = h3;

  if (is_boundary0) {
    h0->setNeighbors(h0->next(), h3, v4, e0, h0->face());
    h10->setNeighbors(h0, h11, v0, e5, h0->face());
  }

  if (is_boundary1) {
    h11->setNeighbors(h3->next(), h10, v4, e5, h3->face());
    h3->setNeighbors(h11, h0, v1, e0, h3->face());
  }

  return v4;
}

VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // This method should collapse the given edge and return an iterator to
  // the new vertex created by the collapse.

  showError("collapseEdge() not implemented.");
  return VertexIter();
}

VertexIter HalfedgeMesh::collapseFace(FaceIter f) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // This method should collapse the given face and return an iterator to
  // the new vertex created by the collapse.
  showError("collapseFace() not implemented.");
  return VertexIter();
}

FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // This method should replace the given vertex and all its neighboring
  // edges and faces with a single face, returning the new face.

  if (v->isBoundary()) {
    showError("eraseVertex() cannot applied to a boundary vertex.");
    return v->halfedge()->face();
  }

  vector<HalfedgeIter> h_ori_out, h_ori_in;
  vector<VertexIter> v_ori;
  vector<EdgeIter> e_ori;
  vector<FaceIter> f_ori;
  FaceIter f_new;

  // Collect original elements
  HalfedgeIter h = v->halfedge();
  do {
    h_ori_out.push_back(h->next());
    HalfedgeIter h_tmp = h;
    do {
      h_tmp = h_tmp->next()->twin();
    } while (h_tmp->next() != h->twin());
    h_ori_in.push_back(h_tmp);
    v_ori.push_back(h->twin()->vertex());
    e_ori.push_back(h->edge());
    f_ori.push_back(h->face());
    h = h->twin()->next();
  } while (h != v->halfedge());

  // Allocate new elements
  f_new = newFace();

  // Reassign elements
  int n = f_ori.size();
  for (int i = 0; i < n; i++) {
    h_ori_in[i]->next() = h_ori_out[i];
    v_ori[i]->halfedge() = h_ori_out[i];
  }
  f_new->halfedge() = h_ori_out[0];
  h = f_new->halfedge();
  do {
    h->face() = f_new;
    h = h->next();
  } while (h != f_new->halfedge());

  for (int i = 0; i < n; i++) {
    deleteHalfedge(e_ori[i]->halfedge()->twin());
    deleteHalfedge(e_ori[i]->halfedge());
    deleteEdge(e_ori[i]);
    deleteFace(f_ori[i]);
  }
  deleteVertex(v);
  return f_new;
}

FaceIter HalfedgeMesh::eraseEdge(EdgeIter e) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // This method should erase the given edge and return an iterator to the
  // merged face.

  if (e->isBoundary()) {
    showError("eraseEdge() cannot applied to boundary edge.");
    if (!e->halfedge()->isBoundary()) {
      return e->halfedge()->face();
    } else {
      return e->halfedge()->twin()->face();
    }
  }

  HalfedgeIter h1, h2, h3, h4;
  VertexIter v1, v2;
  FaceIter f1, f2, f_new;

  // handle edge cases
  v1 = e->halfedge()->twin()->vertex();
  v2 = e->halfedge()->vertex();
  if (v1->degree() <= 1 && v2->degree() <= 1) {
    f1 = e->halfedge()->face();
    deleteHalfedge(e->halfedge()->twin());
    deleteHalfedge(e->halfedge());
    deleteVertex(v1);
    deleteVertex(v2);
    deleteEdge(e);
    return f1;  // this could have error in extreme edge case
  }
  if (v1->degree() <= 1 || v2->degree() <= 1) {
    VertexIter v_tmp;
    HalfedgeIter h_tmp;
    h1 = e->halfedge();
    h2 = e->halfedge()->twin();
    if (v2->degree() <= 1) {
      v_tmp = v2;
      v2 = v1;
      v1 = v_tmp;
      h_tmp = h2;
      h2 = h1;
      h1 = h_tmp;
    }
    h3 = h2->next();
    h_tmp = h2;
    do {
      h_tmp = h_tmp->next()->twin();
    } while (h_tmp->next() != h1);
    h4 = h_tmp;
    h4->next() = h3;
    h3->vertex()->halfedge() = h3;
    h3->face()->halfedge() = h3;
    deleteHalfedge(h1);
    deleteHalfedge(h2);
    deleteVertex(v1);
    deleteEdge(e);
    return h3->face();
  }


  // Collect original elements
  HalfedgeIter h = e->halfedge();
  h1 = h->next();
  do {
    h = h->next()->twin();
  } while (h->next() != e->halfedge()->twin());
  h2 = h;
  h = e->halfedge()->twin();
  h3 = h->next();
  do {
    h = h->next()->twin();
  } while (h->next() != e->halfedge());
  h4 = h;
  v1 = h1->vertex();
  v2 = h3->vertex();
  f1 = e->halfedge()->face();
  f2 = e->halfedge()->twin()->face();

  // Allocate new elements
  f_new = newFace();

  // Reassign elements
  h2->next() = h1;
  h4->next() = h3;
  v1->halfedge() = h1;
  v2->halfedge() = h3;
  f_new->halfedge() = h1;
  h = h1;
  do {
    h->face() = f_new;
    h = h->next();
  } while (h != f_new->halfedge());

  deleteHalfedge(e->halfedge()->twin());
  deleteHalfedge(e->halfedge());
  deleteEdge(e);
  deleteFace(f1);
  deleteFace(f2);
  return f_new;
}

EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
  // TODO: (meshEdit)
  // This method should flip the given edge and return an iterator to the
  // flipped edge.

  // Check boundary
  if (e0->isBoundary()) {
    showError("Cannot flip a boundary edge.");
    return e0;
  }

  HalfedgeIter h, h0, h1, h2, h3, h4, h5, h6, h7;
  VertexIter v0, v1, v2, v3;
  FaceIter f0, f1;

  // Collect original elements
  h0 = e0->halfedge();
  h1 = h0->twin();
  h2 = h0->next();
  h4 = h1->next();
  h6 = h2->next();
  h7 = h4->next();
  v0 = h0->vertex();
  v1 = h1->vertex();
  v2 = h6->vertex();
  v3 = h7->vertex();
  f0 = h0->face();
  f1 = h1->face();

  h = h1;
  do {
    if (h->next() == h1) {
      h3 = h;  // Collect h3
      break;
    }
    h = h->next();
  } while (h != h1);

  h = h0;
  do {
    if (h->next() == h0) {
      h5 = h;  // Collect h5
      break;
    }
    h = h->next();
  } while (h != h0);

  // Reassign elements
  h0->setNeighbors(h6, h1, v3, e0, f0);
  h1->setNeighbors(h7, h0, v2, e0, f1);
  h2->setNeighbors(h1, h2->twin(), h2->vertex(), h2->edge(), f1);
  h3->setNeighbors(h2, h3->twin(), h3->vertex(), h3->edge(), h3->face());
  h4->setNeighbors(h0, h4->twin(), h4->vertex(), h4->edge(), f0);
  h5->setNeighbors(h4, h5->twin(), h5->vertex(), h5->edge(), h5->face());

  v0->halfedge() = h4;
  v1->halfedge() = h2;

  f0->halfedge() = h0;
  f1->halfedge() = h1;

  return e0;
}

void HalfedgeMesh::subdivideQuad(bool useCatmullClark) {
  // Unlike the local mesh operations (like bevel or edge flip), we will perform
  // subdivision by splitting *all* faces into quads "simultaneously."  Rather
  // than operating directly on the halfedge data structure (which as you've
  // seen
  // is quite difficult to maintain!) we are going to do something a bit nicer:
  //
  //    1. Create a raw list of vertex positions and faces (rather than a full-
  //       blown halfedge mesh).
  //
  //    2. Build a new halfedge mesh from these lists, replacing the old one.
  //
  // Sometimes rebuilding a data structure from scratch is simpler (and even
  // more
  // efficient) than incrementally modifying the existing one.  These steps are
  // detailed below.

  // TODO Step I: Compute the vertex positions for the subdivided mesh.  Here
  // we're
  // going to do something a little bit strange: since we will have one vertex
  // in
  // the subdivided mesh for each vertex, edge, and face in the original mesh,
  // we
  // can nicely store the new vertex *positions* as attributes on vertices,
  // edges,
  // and faces of the original mesh.  These positions can then be conveniently
  // copied into the new, subdivided mesh.
  // [See subroutines for actual "TODO"s]
  if (useCatmullClark) {
    computeCatmullClarkPositions();
  } else {
    computeLinearSubdivisionPositions();
  }

  // TODO Step II: Assign a unique index (starting at 0) to each vertex, edge,
  // and
  // face in the original mesh.  These indices will be the indices of the
  // vertices
  // in the new (subdivided mesh).  They do not have to be assigned in any
  // particular
  // order, so long as no index is shared by more than one mesh element, and the
  // total number of indices is equal to V+E+F, i.e., the total number of
  // vertices
  // plus edges plus faces in the original mesh.  Basically we just need a
  // one-to-one
  // mapping between original mesh elements and subdivided mesh vertices.
  // [See subroutine for actual "TODO"s]
  assignSubdivisionIndices();

  // TODO Step III: Build a list of quads in the new (subdivided) mesh, as
  // tuples of
  // the element indices defined above.  In other words, each new quad should be
  // of
  // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
  // original mesh elements.  Note that it is essential to get the orientation
  // right
  // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
  // circulate in the same direction as old faces (think about the right-hand
  // rule).
  // [See subroutines for actual "TODO"s]
  vector<vector<Index> > subDFaces;
  vector<Vector3D> subDVertices;
  buildSubdivisionFaceList(subDFaces);
  buildSubdivisionVertexList(subDVertices);

  // TODO Step IV: Pass the list of vertices and quads to a routine that clears
  // the
  // internal data for this halfedge mesh, and builds new halfedge data from
  // scratch,
  // using the two lists.
  rebuild(subDFaces, subDVertices);
}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * simple linear interpolation, e.g., the edge midpoints and face
 * centroids.
 */
void HalfedgeMesh::computeLinearSubdivisionPositions() {
  // TODO For each vertex, assign Vertex::newPosition to
  // its original position, Vertex::position.
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    v->newPosition = v->position;
  }

  // TODO For each edge, assign the midpoint of the two original
  // positions to Edge::newPosition.
  for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
    VertexIter v1 = e->halfedge()->vertex();
    VertexIter v2 = e->halfedge()->twin()->vertex();
    e->newPosition = (v1->position + v2->position) / 2.0;
  }

  // TODO For each face, assign the centroid (i.e., arithmetic mean)
  // of the original vertex positions to Face::newPosition.  Note
  // that in general, NOT all faces will be triangles!
  for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
    HalfedgeIter h0 = f->halfedge();
    HalfedgeIter h = h0;
    Vector3D sum_point;
    double count = 0.0;
    do {
      count += 1.0;
      sum_point = sum_point + h->vertex()->position;
      h = h->next();
    } while (h != h0);
    f->newPosition = sum_point / count;
  }
}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * the Catmull-Clark rules for subdivision.
 */
void HalfedgeMesh::computeCatmullClarkPositions() {
  // TODO The implementation for this routine should be
  // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
  // except that the calculation of the positions themsevles is
  // slightly more involved, using the Catmull-Clark subdivision
  // rules. (These rules are outlined in the Developer Manual.)

  // TODO face
  for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
    HalfedgeIter h0 = f->halfedge();
    HalfedgeIter h = h0;
    Vector3D sum_point;
    double count = 0.0;
    do {
      count += 1.0;
      sum_point = sum_point + h->vertex()->position;
      h = h->next();
    } while (h != h0);
    f->newPosition = sum_point / count;
  }

  // TODO edges
  for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
    VertexIter v1 = e->halfedge()->vertex();
    VertexIter v2 = e->halfedge()->twin()->vertex();
    if (e->isBoundary()) {
      e->newPosition = (v1->position + v2->position) / 2.0;
    } else {
      FaceIter f1 = e->halfedge()->face();
      FaceIter f2 = e->halfedge()->twin()->face();
      e->newPosition = (v1->position + v2->position + f1->newPosition + f2->newPosition) / 4.0;
    }
  }

  // TODO vertices
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    HalfedgeIter h0 = v->halfedge();
    HalfedgeIter h = h0;
    if (v->isBoundary()) {
      double count = 0.0;
      Vector3D sum_point;
      do {
        EdgeIter e = h->edge();
        if (e->isBoundary()) {
          count += 1.0;
          sum_point = sum_point + e->newPosition;
        }
        h = h->twin()->next();
      } while (h != h0);
      if (count > 0.0) {
        v->newPosition = sum_point / (count * 4.0) + v->position * 3.0 / 4.0;
      } else {
        v->newPosition = v->position;
      }
    } else {
      double n = 0.0;
      Vector3D Q, R, S;
      do {
        FaceIter f = h->face();
        EdgeIter e = h->edge();
        n += 1.0;
        Q = Q + f->newPosition;
        R = R + e->newPosition;
        h = h->twin()->next();
      } while (h != h0);
      Q = Q / n;
      R = R / n;
      S = v->position;
      if (n >= 3.0) {
        v->newPosition = (Q + R * 2.0 + S * (n - 3.0)) / n;
      } else {
        v->newPosition = (Q + R * 2.0) / 3.0;
      }
    }
  }
}

/**
 * Assign a unique integer index to each vertex, edge, and face in
 * the mesh, starting at 0 and incrementing by 1 for each element.
 * These indices will be used as the vertex indices for a mesh
 * subdivided using Catmull-Clark (or linear) subdivision.
 */
void HalfedgeMesh::assignSubdivisionIndices() {
  // TODO Start a counter at zero; if you like, you can use the
  // "Index" type (defined in halfedgeMesh.h)
  Index counter = 0;

  // TODO Iterate over vertices, assigning values to Vertex::index
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    v->index = counter;
    counter++;
  }

  // TODO Iterate over edges, assigning values to Edge::index
  for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
    e->index = counter;
    counter++;
  }

  // TODO Iterate over faces, assigning values to Face::index
  for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
    f->index = counter;
    counter++;
  }
}

/**
 * Build a flat list containing all the vertex positions for a
 * Catmull-Clark (or linear) subdivison of this mesh.  The order of
 * vertex positions in this list must be identical to the order
 * of indices assigned to Vertex::newPosition, Edge::newPosition,
 * and Face::newPosition.
 */
void HalfedgeMesh::buildSubdivisionVertexList(vector<Vector3D>& subDVertices) {
  // TODO Resize the vertex list so that it can hold all the vertices.
  subDVertices.reserve(nVertices() + nEdges() + nFaces());

  // TODO Iterate over vertices, assigning Vertex::newPosition to the
  // appropriate location in the new vertex list.
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    subDVertices.push_back(v->newPosition);
  }

  // TODO Iterate over edges, assigning Edge::newPosition to the appropriate
  // location in the new vertex list.
  for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
    subDVertices.push_back(e->newPosition);
  }

  // TODO Iterate over faces, assigning Face::newPosition to the appropriate
  // location in the new vertex list.
  for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
    subDVertices.push_back(f->newPosition);
  }
}

/**
 * Build a flat list containing all the quads in a Catmull-Clark
 * (or linear) subdivision of this mesh.  Each quad is specified
 * by a vector of four indices (i,j,k,l), which come from the
 * members Vertex::index, Edge::index, and Face::index.  Note that
 * the ordering of these indices is important because it determines
 * the orientation of the new quads; it is also important to avoid
 * "bowties."  For instance, (l,k,j,i) has the opposite orientation
 * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
 * will look like a bowtie.
 */
void HalfedgeMesh::buildSubdivisionFaceList(vector<vector<Index> >& subDFaces) {
  // TODO This routine is perhaps the most tricky step in the construction of
  // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
  // vertex positions).  Basically what you want to do is iterate over faces,
  // then for each for each face, append N quads to the list (where N is the
  // degree of the face).  For this routine, it may be more convenient to simply
  // append quads to the end of the list (rather than allocating it ahead of
  // time), though YMMV.  You can of course iterate around a face by starting
  // with its first halfedge and following the "next" pointer until you get
  // back to the beginning.  The tricky part is making sure you grab the right
  // indices in the right order---remember that there are indices on vertices,
  // edges, AND faces of the original mesh.  All of these should get used.  Also
  // remember that you must have FOUR indices per face, since you are making a
  // QUAD mesh!

  // TODO iterate over faces
  // TODO loop around face
  // TODO build lists of four indices for each sub-quad
  // TODO append each list of four indices to face list
  for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
    HalfedgeIter h0 = f->halfedge();
    HalfedgeIter h, h_next;
    EdgeIter e, e_next;
    VertexIter v;

    h = h0;
    e = h->edge();
    do {
      h_next = h->next();
      e_next = h_next->edge();
      v = h_next->vertex();

      vector<Index> quad(4);
      quad[0] = e->index;
      quad[1] = v->index;
      quad[2] = e_next->index;
      quad[3] = f->index;
      subDFaces.push_back(quad);

      h = h_next;
      e = e_next;
    } while (h != h0);
  }
}

FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {
  // *** Extra Credit ***
  // TODO This method should replace the vertex v with a face, corresponding to
  // a bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelVertexComputeNewPositions (which you also have to
  // implement!)

  if (v->isBoundary()) {
    showError("bevelVertex() cannot applied to boundary vertex.");
    return facesBegin();
  }

  if (v->degree() < 3) {
    showError("bevelVertex() cannot applied to vertex of degree less than 3.");
    return facesBegin();
  }

  int n = v->degree();
  vector<HalfedgeIter> h_ori_out, h_ori_in;
  vector<FaceIter> f_ori;
  vector<HalfedgeIter> h_new_inter, h_new_outer;
  vector<VertexIter> v_new;
  vector<EdgeIter> e_new;
  FaceIter f_new;

  // Collect original elements
  HalfedgeIter h = v->halfedge();
  do {
    h_ori_out.push_back(h);
    h_ori_in.push_back(h->twin());
    f_ori.push_back(h->twin()->face());
    h = h->twin()->next();
  } while (h != v->halfedge());

  if (f_ori.size() != n) {
    showError("vertex degree not equal to associated faces number.");
    return facesBegin();
  }

  // Allocate new elements
  for (int i = 0; i < n; i++) {
    h_new_inter.push_back(newHalfedge());
    h_new_outer.push_back(newHalfedge());
    v_new.push_back(newVertex());
    e_new.push_back(newEdge());
  }
  f_new = newFace();

  // Reassign elements
  for (int i = 0; i < n; i++) {
    h_ori_out[i]->vertex() = v_new[i];
    h_ori_in[i]->next() = h_new_outer[i];
    h_new_inter[i]->setNeighbors(h_new_inter[(i - 1 + n) % n], h_new_outer[i],
      v_new[(i + 1) % n], e_new[i], f_new);
    h_new_outer[i]->setNeighbors(h_ori_out[(i + 1) % n], h_new_inter[i],
      v_new[i], e_new[i], f_ori[i]);
    v_new[i]->halfedge() = h_new_outer[i];
    e_new[i]->halfedge() = h_new_outer[i];

    v_new[i]->position = v->position;
  }
  f_new->halfedge() = h_new_inter[0];

  deleteVertex(v);
  return f_new;
}

FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {
  // *** Extra Credit ***
  // TODO This method should replace the edge e with a face, corresponding to a
  // bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelEdgeComputeNewPositions (which you also have to
  // implement!)

  if (e->halfedge()->vertex()->isBoundary() || e->halfedge()->twin()->vertex()->isBoundary()) {
    showError("bevelEdge() cannot applied to edge with boundary vertex.");
    return facesBegin();
  }

  int n = 0, m;
  vector<HalfedgeIter> h_ori_out, h_ori_in;
  vector<FaceIter> f_ori;
  vector<HalfedgeIter> h_new_inter, h_new_outer;
  vector<VertexIter> v_new;
  vector<EdgeIter> e_new;
  FaceIter f_new;

  // Collect original elements
  HalfedgeIter h = e->halfedge()->next();
  do {
    n++;
    h_ori_out.push_back(h);
    h_ori_in.push_back(h->twin());
    f_ori.push_back(h->twin()->face());
    h = h->twin()->next();
  } while (h != e->halfedge()->twin());
  m = n;
  h = h->next();
  do {
    n++;
    h_ori_out.push_back(h);
    h_ori_in.push_back(h->twin());
    f_ori.push_back(h->twin()->face());
    h = h->twin()->next();
  } while (h != e->halfedge());

  // Allocate new elements
  for (int i = 0; i < n; i++) {
    h_new_inter.push_back(newHalfedge());
    h_new_outer.push_back(newHalfedge());
    v_new.push_back(newVertex());
    e_new.push_back(newEdge());
  }
  f_new = newFace();

  // Reassign elements
  VertexIter v1 = e->halfedge()->twin()->vertex();
  VertexIter v2 = e->halfedge()->vertex();
  for (int i = 0; i < n; i++) {
    h_ori_out[i]->vertex() = v_new[i];
    h_ori_in[i]->next() = h_new_outer[i];
    f_ori[i]->halfedge() = h_ori_in[i];
    h_new_inter[i]->setNeighbors(h_new_inter[(i - 1 + n) % n], h_new_outer[i],
      v_new[(i + 1) % n], e_new[i], f_new);
    h_new_outer[i]->setNeighbors(h_ori_out[(i + 1) % n], h_new_inter[i],
      v_new[i], e_new[i], f_ori[i]);
    v_new[i]->halfedge() = h_new_outer[i];
    e_new[i]->halfedge() = h_new_outer[i];

    if (i < m) {
      v_new[i]->position = v1->position;
    } else {
      v_new[i]->position = v2->position;
    }
  }
  f_new->halfedge() = h_new_inter[0];

  deleteVertex(v1);
  deleteVertex(v2);
  deleteEdge(e);
  return f_new;
}

FaceIter HalfedgeMesh::bevelFace(FaceIter f) {
  // *** Extra Credit ***
  // TODO This method should replace the face f with an additional, inset face
  // (and ring of faces around it), corresponding to a bevel operation. It
  // should return the new face.  NOTE: This method is responsible for updating
  // the *connectivity* of the mesh only---it does not need to update the vertex
  // positions.  These positions will be updated in
  // HalfedgeMesh::bevelFaceComputeNewPositions (which you also have to
  // implement!)

  if (f->isBoundary()) {
    showError("bevelFace() cannot applied to boundary face.");
    return facesBegin();
  }

  if (f->degree() < 3) {
    showError("bevelFace() cannot applied to face of degree less than 3.");
    return facesBegin();
  }

  int n = f->degree();
  vector<HalfedgeIter> h_ori;
  vector<VertexIter> v_ori;
  FaceIter f_new_center;
  vector<HalfedgeIter> h_new_inter, h_new_outer, h_new_incoming, h_new_outgoing;
  vector<VertexIter> v_new;
  vector<EdgeIter> e_new, e_new_bridge;
  vector<FaceIter> f_new;
  Vector3D center_p;

  // Collect original elements
  HalfedgeIter h = f->halfedge();
  do {
    h_ori.push_back(h);
    v_ori.push_back(h->vertex());
    center_p = center_p + h->vertex()->position;
    h = h->next();
  } while (h != f->halfedge());
  center_p = center_p / double(n);

  if (v_ori.size() != n) {
    showError("face degree not equal to associated vertices number.");
    return facesBegin();
  }

  // Allocate new elements
  for (int i = 0; i < n; i++) {
    h_new_inter.push_back(newHalfedge());
    h_new_outer.push_back(newHalfedge());
    h_new_incoming.push_back(newHalfedge());
    h_new_outgoing.push_back(newHalfedge());
    v_new.push_back(newVertex());
    e_new.push_back(newEdge());
    e_new_bridge.push_back(newEdge());
    f_new.push_back(newFace());
  }
  f_new_center = newFace();

  // Reassign elements
  for (int i = 0; i < n; i++) {
    h_ori[i]->next() = h_new_incoming[(i + 1) % n];
    h_ori[i]->face() = f_new[i];
    h_new_inter[i]->setNeighbors(h_new_inter[(i + 1) % n], h_new_outer[i],
      v_new[i], e_new[i], f_new_center);
    h_new_outer[i]->setNeighbors(h_new_outgoing[i], h_new_inter[i],
      v_new[(i + 1) % n], e_new[i], f_new[i]);
    h_new_incoming[i]->setNeighbors(h_new_outer[(i - 1 + n) % n], h_new_outgoing[i],
      v_ori[i], e_new_bridge[i], f_new[(i - 1 + n) % n]);
    h_new_outgoing[i]->setNeighbors(h_ori[i], h_new_incoming[i],
      v_new[i], e_new_bridge[i], f_new[i]);
    v_new[i]->halfedge() = h_new_inter[i];
    e_new[i]->halfedge() = h_new_inter[i];
    e_new_bridge[i]->halfedge() = h_new_outgoing[i];
    f_new[i]->halfedge() = h_new_outgoing[i];

    v_new[i]->position = (v_ori[i]->position + center_p) / 2.0;
  }
  f_new_center->halfedge() = h_new_inter[0];

  deleteFace(f);
  return f_new_center;
}


void HalfedgeMesh::bevelFaceComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double normalShift,
    double tangentialInset) {
  // *** Extra Credit ***
  // TODO Compute new vertex positions for the vertices of the beveled face.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., newHalfedges.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the originalVertexPositions array) to compute an offset vertex
  // position.
  //
  // Note that there is a 1-to-1 correspondence between halfedges in
  // newHalfedges and vertex positions
  // in orig.  So, you can write loops of the form
  //
  // for( int i = 0; i < newHalfedges.size(); hs++ )
  // {
  //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
  //    position correponding to vertex i
  // }
  //

  // Compute tangential offset ratio t
  Vector3D ori_center_p;
  for (auto ori_p : originalVertexPositions) {
    ori_center_p = ori_center_p + ori_p;
  }
  ori_center_p = ori_center_p / double(originalVertexPositions.size());

  Vector3D new_center_p;
  for (auto h : newHalfedges) {
    new_center_p = new_center_p + h->vertex()->position;
  }
  new_center_p = new_center_p / double(newHalfedges.size());

  Vector3D v_ori_p = originalVertexPositions[0];
  Vector3D v_new_p = newHalfedges[0]->vertex()->position;
  double t = (v_new_p - new_center_p).norm() / (v_ori_p - ori_center_p).norm();
  t = t + tangentialInset * 2;
  t = std::max(t, 1e-6);

  // Compute normal offset ratio r
  Vector3D normal_vec = newHalfedges[0]->twin()->next()->twin()->face()->normal();
  Vector3D normal_offset = new_center_p - ori_center_p;
  if (normal_vec.x == 0.0 && normal_vec.y == 0.0 && normal_vec.z == 0.0) {
    showError("error: normal vector of face is (0, 0, 0).");
    return;
  }
  double r;
  if (std::abs(normal_vec.x) > 1e-1) {
    r = normal_offset.x / normal_vec.x;
  } else if (std::abs(normal_vec.y) > 1e-1) {
    r = normal_offset.y / normal_vec.y;
  } else {
    r = normal_offset.z / normal_vec.z;
  }
  r = r - normalShift * 10;

  // Reassign new vertex positions
  for (int i = 0; i < newHalfedges.size(); i++) {
    newHalfedges[i]->vertex()->position = ori_center_p + (originalVertexPositions[i] - ori_center_p) * t + normal_vec * r;
  }
}

void HalfedgeMesh::bevelVertexComputeNewPositions(
    Vector3D originalVertexPosition, vector<HalfedgeIter>& newHalfedges,
    double tangentialInset) {
  // *** Extra Credit ***
  // TODO Compute new vertex positions for the vertices of the beveled vertex.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., hs.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.

  VertexIter v_new = newHalfedges[0]->vertex();
  VertexIter v_ori = newHalfedges[0]->twin()->vertex();
  Vector3D center_p = originalVertexPosition;
  double t = (v_new->position - center_p).norm() / (v_ori->position - center_p).norm();
  t = t + tangentialInset * 2;
  t = std::max(std::min(t, 1.0), 0.0);

  int n = newHalfedges.size();
  for (int i = 0; i < n; i++) {
    v_new = newHalfedges[i]->vertex();
    v_ori = newHalfedges[i]->twin()->vertex();
    v_new->position = center_p + (v_ori->position - center_p) * t;
  }
}

void HalfedgeMesh::bevelEdgeComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double tangentialInset) {
  // *** Extra Credit ***
  // TODO Compute new vertex positions for the vertices of the beveled edge.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., newHalfedges.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.
  //
  // Note that there is a 1-to-1 correspondence between halfedges in
  // newHalfedges and vertex positions
  // in orig.  So, you can write loops of the form
  //
  // for( int i = 0; i < newHalfedges.size(); i++ )
  // {
  //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
  //    position correponding to vertex i
  // }
  //

  VertexIter v_new = newHalfedges[0]->vertex();
  VertexIter v_ori = newHalfedges[0]->twin()->vertex();
  Vector3D center_p = originalVertexPositions[0];
  double t = (v_new->position - center_p).norm() / (v_ori->position - center_p).norm();
  t = t + tangentialInset * 2;
  t = std::max(std::min(t, 1.0), 0.0);

  int n = newHalfedges.size();
  for (int i = 0; i < n; i++) {
    v_new = newHalfedges[i]->vertex();
    v_ori = newHalfedges[i]->twin()->vertex();
    center_p = originalVertexPositions[i];
    v_new->position = center_p + (v_ori->position - center_p) * t;
  }
}

void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
  for (auto f : fcs) splitPolygon(f);
}

void HalfedgeMesh::splitPolygon(FaceIter f) {
  // *** Extra Credit ***
  // TODO: (meshedit) 
  // Triangulate a polygonal face
  int count = 0;
  HalfedgeIter h = f->halfedge();
  do {
    count++;
    h = h->next();
  } while (h != f->halfedge());
  if (count <= 3) {
    return;
  }

  FaceIter f0, f1;
  HalfedgeIter h0, h1, h2, h3, h4, h5;
  VertexIter v0, v1;
  EdgeIter e0;
  
  // Collect original elements
  f0 = f;
  h0 = f0->halfedge();
  h1 = h0->next();
  h2 = h1->next();
  h = h0;
  do {
    h = h->next();
  } while (h->next() != h0);
  h3 = h;
  v0 = h0->vertex();
  v1 = h2->vertex();

  // Allocate new elements
  h4 = newHalfedge();
  h5 = newHalfedge();
  e0 = newEdge();
  f1 = newFace();

  // Reassign elements
  h1->next() = h4;
  h2->face() = f1;
  h3->next() = h5;
  h3->face() = f1;
  h4->setNeighbors(h0, h5, v1, e0, f0);
  h5->setNeighbors(h2, h4, v0, e0, f1);
  e0->halfedge() = h4;
  f1->halfedge() = h5;

  splitPolygon(f1);
}

EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // Compute the combined quadric from the edge endpoints.
  // -> Build the 3x3 linear system whose solution minimizes the quadric error
  //    associated with these two endpoints.
  // -> Use this system to solve for the optimal position, and store it in
  //    EdgeRecord::optimalPoint.
  // -> Also store the cost associated with collapsing this edg in
  //    EdgeRecord::Cost.
}

void MeshResampler::upsample(HalfedgeMesh& mesh)
// This routine should increase the number of triangles in the mesh using Loop
// subdivision.
{
  // TODO: (meshEdit)
  // Compute new positions for all the vertices in the input mesh, using
  // the Loop subdivision rule, and store them in Vertex::newPosition.
  // -> At this point, we also want to mark each vertex as being a vertex of the
  //    original mesh.
  // -> Next, compute the updated vertex positions associated with edges, and
  //    store it in Edge::newPosition.
  // -> Next, we're going to split every edge in the mesh, in any order.  For
  //    future reference, we're also going to store some information about which
  //    subdivided edges come from splitting an edge in the original mesh, and
  //    which edges are new, by setting the flat Edge::isNew. Note that in this
  //    loop, we only want to iterate over edges of the original mesh.
  //    Otherwise, we'll end up splitting edges that we just split (and the
  //    loop will never end!)
  // -> Now flip any new edge that connects an old and new vertex.
  // -> Finally, copy the new vertex positions into final Vertex::position.

  // Each vertex and edge of the original surface can be associated with a
  // vertex in the new (subdivided) surface.
  // Therefore, our strategy for computing the subdivided vertex locations is to
  // *first* compute the new positions
  // using the connectity of the original (coarse) mesh; navigating this mesh
  // will be much easier than navigating
  // the new subdivided (fine) mesh, which has more elements to traverse.  We
  // will then assign vertex positions in
  // the new mesh based on the values we computed for the original mesh.

  // Check triangle meshes
  for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
    HalfedgeIter h0 = f->halfedge();
    HalfedgeIter h = h0;
    int count = 0;
    do {
      count++;
      h = h->next();
    } while (h != h0);
    if (count != 3) {
      showError("upsample() cannot be applied to non-triangle meshes.");
      return;
    }
  }

  // Mark all original vertices as old
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    v->isNew = false;
  }

  // Next, compute the updated vertex positions associated with edges.
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    VertexIter v1 = e->halfedge()->vertex();
    VertexIter v2 = e->halfedge()->twin()->vertex();
    if (e->isBoundary()) {
      e->newPosition = (v1->position + v2->position) / 2.0;
    } else {
      VertexIter v3 = e->halfedge()->next()->next()->vertex();
      VertexIter v4 = e->halfedge()->twin()->next()->next()->vertex();
      e->newPosition = (v1->position + v2->position) * 3.0 / 8.0 + (v3->position + v4->position) / 8.0;
    }
  }

  // Compute updated positions for all the vertices in the original mesh, using
  // the Loop subdivision rule.
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    HalfedgeIter h0 = v->halfedge();
    HalfedgeIter h = h0;
    if (v->isBoundary()) {
      double count = 0.0;
      Vector3D sum_point;
      do {
        EdgeIter e = h->edge();
        if (e->isBoundary()) {
          count += 1.0;
          sum_point = sum_point + e->newPosition;
        }
        h = h->twin()->next();
      } while (h != h0);
      if (count > 0.0) {
        v->newPosition = sum_point / (count * 4.0) + v->position * 3.0 / 4.0;
      } else {
        v->newPosition = v->position;
      }
    } else {
      double n = 0.0;
      Vector3D sum_point;
      do {
        VertexIter v = h->twin()->vertex();
        n += 1.0;
        sum_point = sum_point + v->position;
        h = h->twin()->next();
      } while (h != h0);
      if (n > 3.0) {
        v->newPosition = v->position * 5.0 / 8.0 + sum_point * 3.0 / (8.0 * n);
      } else if (n > 0.0) {
        v->newPosition = v->position * (16.0 - n * 3.0) / 16.0 + sum_point * 3.0 / 16.0;
      } else {
        v->newPosition = v->position;
      }
    }
  }

  // Next, we're going to split every edge in the mesh, in any order.  For
  // future
  // reference, we're also going to store some information about which
  // subdivided
  // edges come from splitting an edge in the original mesh, and which edges are
  // new.
  // In this loop, we only want to iterate over edges of the original
  // mesh---otherwise,
  // we'll end up splitting edges that we just split (and the loop will never
  // end!)
  int num_edges = mesh.nEdges();
  EdgeIter e = mesh.edgesBegin();
  for (int i = 0; i < num_edges; i++) {
    EdgeIter nextEdge = e;
    nextEdge++;

    VertexIter v1 = e->halfedge()->vertex();
    VertexIter v2 = e->halfedge()->twin()->vertex();
    Vector3D e_newPosition = e->newPosition * 1.0;

    // Check for edge cases
    if (e->halfedge()->isBoundary() && e->halfedge()->twin()->isBoundary()) {
      e->isNew = false;
      e = nextEdge;
      continue;
    }

    // Split edge and set information
    VertexIter v_new = mesh.splitEdge(e);
    v_new->isNew = true;
    v_new->newPosition = e_newPosition;

    HalfedgeIter h0 = v_new->halfedge();
    HalfedgeIter h = h0;
    do {
      VertexIter v3 = h->twin()->vertex();
      EdgeIter e3 = h->edge();
      if (v3 == v1 || v3 == v2) {
        e3->isNew = false;
      } else {
        e3->isNew = true;
      }
      h = h->twin()->next();
    } while (h != h0);

    e = nextEdge;
  }

  // Finally, flip any new edge that connects an old and new vertex.
  for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
    if (e->isBoundary()) {
      continue;
    }

    VertexIter v1 = e->halfedge()->vertex();
    VertexIter v2 = e->halfedge()->twin()->vertex();
    if (e->isNew && (v1->isNew ^ v2->isNew)) {
      mesh.flipEdge(e);
    }
  }

  // Copy the updated vertex positions to the subdivided mesh.
  for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
    v->position = v->newPosition;
  }
}

void MeshResampler::downsample(HalfedgeMesh& mesh) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // Compute initial quadrics for each face by simply writing the plane equation
  // for the face in homogeneous coordinates. These quadrics should be stored
  // in Face::quadric
  // -> Compute an initial quadric for each vertex as the sum of the quadrics
  //    associated with the incident faces, storing it in Vertex::quadric
  // -> Build a priority queue of edges according to their quadric error cost,
  //    i.e., by building an EdgeRecord for each edge and sticking it in the
  //    queue.
  // -> Until we reach the target edge budget, collapse the best edge. Remember
  //    to remove from the queue any edge that touches the collapsing edge
  //    BEFORE it gets collapsed, and add back into the queue any edge touching
  //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
  //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
  //    top of the queue.
  showError("downsample() not implemented.");
}

void MeshResampler::resample(HalfedgeMesh& mesh) {
  // *** Extra Credit ***
  // TODO: (meshEdit)
  // Compute the mean edge length.
  // Repeat the four main steps for 5 or 6 iterations
  // -> Split edges much longer than the target length (being careful about
  //    how the loop is written!)
  // -> Collapse edges much shorter than the target length.  Here we need to
  //    be EXTRA careful about advancing the loop, because many edges may have
  //    been destroyed by a collapse (which ones?)
  // -> Now flip each edge if it improves vertex degree
  // -> Finally, apply some tangential smoothing to the vertex positions
  showError("resample() not implemented.");
}

}  // namespace CS248
