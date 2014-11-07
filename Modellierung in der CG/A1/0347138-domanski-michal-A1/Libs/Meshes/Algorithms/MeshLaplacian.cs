
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MathNet.Numerics.Algorithms.LinearAlgebra;
using Meshes.Generic;
using Meshes.Utilities;
using SharpDX;

using CSparse.Double;

using System.Runtime.Serialization.Formatters.Binary;
using System.IO;

namespace Meshes.Algorithms
{
    /// <summary>
    /// This class encapsulates static methods needed for Laplacian mesh processing. 
    /// </summary>
    public static class MeshLaplacian
    {
        private const double EPSILON = 0.0000000001d;

        public enum Type
        {
            Uniform, Cotan, MeanValue,
        }

        /// <summary>
        /// Collection of supported Laplacina types
        /// </summary>
        public static IList<string> TypeCollection = new List<string>()
        {
            Type.Uniform.ToString(),
            Type.Cotan.ToString(),
            Type.MeanValue.ToString(),                
        };

        /// <summary>
        /// Returns a Laplacian Matrix depending on the chosen type.
        /// </summary>
        public static TripletMatrix CreateLaplacian(TriangleMesh mesh, double lambda = 0.0, double eye = 0.0)
        {
            switch (SelectedLaplacian)
            {
                default:
                case Type.Uniform:
                    return CreateUniformLaplacian(mesh, lambda, eye, LaplacianNormalize);
                case Type.Cotan:
                    return CreateCotanLaplacian(mesh, lambda, eye, LaplacianNormalize);
                case Type.MeanValue:
                    return CreateMeanValueLaplacian(mesh, lambda, eye, LaplacianNormalize);
            }            
        }

        /// <summary>
        /// Create an extended Laplacina with a additinoal constraint matrix appended on the bottom.         
        /// </summary>
        public static TripletMatrix CreateExtendedLaplacian(TriangleMesh mesh, double lambda = 1.0)
        {
            switch (SelectedLaplacian)
            {
                default:
                case Type.MeanValue:
                case Type.Uniform:
                    return CreateExtendedUniformLaplacian(mesh, lambda, LaplacianNormalize);
                case Type.Cotan:
                    return CreateExtendedCotanLaplacian(mesh, lambda, LaplacianNormalize);                                                                     
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public static TripletMatrix CreateAdjacenceMatrix(TriangleMesh mesh)
        {
            var n = mesh.Vertices.Count;
            int nz = mesh.Vertices.Aggregate(0, (c, x) => c += x.VertexCount());
            var A = new TripletMatrix(n, n, nz + n);

            /// build adjacence matrix
            for (int i = 0; i < n; i++)
            {
                var v = mesh.Vertices[i];                
                foreach (var vj in v.Vertices)
                {
                    int j = vj.Index;
                    A[i, j] = 1;                    
                }
            }
            return A;
        }

        /// <summary>
        /// Creates a square uniform n x n laplacian matrix
        /// eye*[ I ] + lambda*[ L ];
        /// Note that depending on lambda you can use this function to create both:
        /// [I - λL]
        /// [I + λL]
        /// If the normalized flag is true, each row of L is normalized to sum to 1
        /// </summary>       
        public static TripletMatrix CreateUniformLaplacian(TriangleMesh mesh, double lambda = 0.0, double eye = 0.0, bool normalized = false)
        {
            var vertexCount = mesh.Vertices.Count;
            int neighborCount = mesh.Vertices.Aggregate(0, (c, x) => x.VertexCount());
            var L = new TripletMatrix(vertexCount, vertexCount, neighborCount, true);

            foreach (var currentVertex in mesh.Vertices)
            {
                // Diagonal entry v(i,i) = eye*1 + λ*L(i,i)
                var diagonalVal = eye + GetNormalized(currentVertex.Vertices.Count(), currentVertex.VertexCount(), normalized) * lambda;
                L.Entry(currentVertex.Index, currentVertex.Index, diagonalVal);
                currentVertex.Vertices.Apply(nb => 
                    L.Entry(currentVertex.Index, nb.Index, GetNormalized(-1, currentVertex.VertexCount(), normalized) * lambda));
            }

            return L;
        }

        private static double GetNormalized(decimal value, decimal sumOfValues, bool normalized)
        {
            return (double) (normalized ? value/sumOfValues : value);
        }

        /// <summary>
        /// Creates a square cotangents n x n laplacian matrix:
        /// eye*[ I ] - lambda*[ L ];
        /// Note that depending on lambda you can use this function to create both:
        /// [I - λL]
        /// [I + λL]
        /// If the normalized flag is true, each row of L is normalized to sum to 1
        /// </summary>        
        public static TripletMatrix CreateCotanLaplacian(TriangleMesh mesh, double lambda = 1.0, double eye = 0.0, bool normalized = false)
        {
            var vertexCount = mesh.Vertices.Count;
            var neighbourCount = mesh.Vertices.Aggregate(0, (c, x) => c += x.VertexCount());
            var L = new TripletMatrix(vertexCount, vertexCount, neighbourCount + vertexCount);

            PrecomputeTraits(mesh);
            foreach (var currentVertex in mesh.Vertices)
            {
                var sumOfAreas = 0d;
                var sumOfWeights = 0d;

                foreach (var halfEdge in currentVertex.Halfedges)
                {
                    sumOfAreas += halfEdge.Traits.VoronoiRegionArea;
                    sumOfWeights += halfEdge.Traits.Cotan;
                }

                var normalization = normalized ? sumOfWeights : sumOfAreas;

                L.Entry(currentVertex.Index, currentVertex.Index, eye + lambda * sumOfWeights / normalization);

                foreach (var halfEdge in currentVertex.Halfedges)
                {
                    L.Entry(currentVertex.Index, halfEdge.ToVertex.Index, lambda * -halfEdge.Traits.Cotan / normalization);
                }
            }
      
            return L;
        }

        /// <summary>
        /// 
        /// This method is not relevant for Assignment 1 yet!
        /// 
        /// Creates a square mean-value n x n laplacian matrix:
        /// eye*[ I ] - lambda*[ L ];
        /// Note that depending on lambda you can use this function to create both:
        /// [I - λL]
        /// [I + λL]
        /// If the normalized flag is true, each row of L is normalized to sum to 1
        /// </summary>
        public static TripletMatrix CreateMeanValueLaplacian(TriangleMesh mesh, double lambda = 1.0, double eye = 0.0, bool normalized = false)
        {
            var n = mesh.Vertices.Count;
            int nz = mesh.Vertices.Aggregate(0, (c, x) => c += x.VertexCount());
            var L = new TripletMatrix(n, n, nz + n);

            return L;
        }

        /// <summary>
        /// Creates an extended uniform n x n laplacian matrix extended by a n x n weight-matrix:
        ///        [ L ]
        /// lambda*[ I ]             
        /// If the normalized flag is true, each row of L is normalized to sum to 1
        /// </summary>       
        public static TripletMatrix CreateExtendedUniformLaplacian(TriangleMesh mesh, double lambda = 1.0, bool normalized = true)
        {
            var vertexCount = mesh.Vertices.Count;

            var L = CreateUniformLaplacian(mesh, 1f, 0f, normalized);
            var identity = new TripletMatrix(vertexCount, vertexCount, vertexCount).SetIdentity();
            identity.MultiplyBy(lambda);

            return L.ConcatenateBelow(identity);
        }

        /// <summary>
        /// Creates a cotangents n x n laplacian matrix extended by a n x n weight-matrix:
        ///        [ L ]
        /// lambda*[ I ]        
        /// If the normalized flag is true, each row of L is normalized to sum to 1
        /// </summary>        
        public static TripletMatrix CreateExtendedCotanLaplacian(TriangleMesh mesh, double lambda = 1.0, bool normalized = true)
        {
            var vertexCount = mesh.Vertices.Count;

            var L = CreateCotanLaplacian(mesh, 1f, 0f, normalized);
            var identity = new TripletMatrix(vertexCount, vertexCount, vertexCount).SetIdentity();
            identity.MultiplyBy(lambda);

            return L.ConcatenateBelow(identity);
        }

        /// <summary>
        /// Extracts the original mesh world-space coordinates as 3 arrays
        /// </summary>
        public static void GetEuclideanCoordinates(TriangleMesh mesh, out double[] px, out double[] py, out double[] pz)
        {
            /// original mesh world-space coordinates
            px = mesh.Vertices.Select(a => (double)a.Traits.Position.X).ToArray();
            py = mesh.Vertices.Select(a => (double)a.Traits.Position.Y).ToArray();
            pz = mesh.Vertices.Select(a => (double)a.Traits.Position.Z).ToArray();
        }

        /// <summary>
        /// Helper function: applies the matrix L to the vectors px, py and pz separately and returns the resutls.
        /// Can compute the differential (=Laplacian) coordinates of the mesh as 3 arrays
        /// </summary>
        public static void ComputeDifferentialCoordinates(SparseMatrix L, double[] px, double[] py, double[] pz,
            out double[] dx, out double[] dy, out double[] dz)
        {
            /// multiply dx = L*px
            L.Ax(px, out dx);
            L.Ax(py, out dy);
            L.Ax(pz, out dz);
        }

        /// <summary>
        /// Update the vertex positions in the mesh
        /// </summary>
        public static void UpdateMesh(TriangleMesh mesh, double[] xx, double[] xy, double[] xz)
        {
            /// update mesh           
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                mesh.Vertices[i].Traits.Position = new Vector3((float)xx[i], (float)xy[i], (float)xz[i]);
            }
        }

        /// <summary>
        /// Update the vertex uv-coordinates in the mesh
        /// </summary>
        public static void UpdateMesh(TriangleMesh mesh, double[] u, double[] v)
        {
            /// update mesh           
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {                
                mesh.Vertices[i].Traits.TextureCoordinate = new Vector2((float)u[i], (float)v[i]);
            }
        }

        /// <summary>
        /// Update the vertex positions and uv-coordinates in the mesh
        /// </summary>
        public static void UpdateMesh(TriangleMesh mesh, double[] xx, double[] xy, double[] xz, double[] u, double[] v)
        {
            /// update mesh           
            for (int i = 0; i < mesh.Vertices.Count; i++)
            {
                mesh.Vertices[i].Traits.Position = new Vector3((float)xx[i], (float)xy[i], (float)xz[i]);
                mesh.Vertices[i].Traits.TextureCoordinate = new Vector2((float)u[i], (float)v[i]);
            }
        }

        /// <summary>
        /// Precompute Traits for all Mesh-Traits that are needed to the cotan-laplacian
        /// Store each angle-cotan at the (outgoing) halfedge of each vertex.
        /// Store each voronoi-area at the (outgoing) halfedge of 
        /// </summary>
        private static void PrecomputeTraits(TriangleMesh mesh)
        {
            /// compute all cotans
            /// store each cotan at the foot of each halfedge
            foreach (var face in mesh.Faces)
            {

                /// compute cotans
                foreach (var halfedge in face.Halfedges)
                {
                    var xiToY = halfedge.ToVertex.Traits.Position - halfedge.FromVertex.Traits.Position;
                    var xiToYprev = halfedge.Opposite.Previous.FromVertex.Traits.Position - halfedge.FromVertex.Traits.Position;
                    var xiToYnext = halfedge.Previous.FromVertex.Traits.Position - halfedge.FromVertex.Traits.Position;
                    var yToYNext = xiToYnext - xiToY;
                    var yPrevToY = xiToY - xiToYprev;

                    // clamp value to 0. As This would correspond to the point where the triangle gets obstruse,
                    // the point effectively gets clamped to the edge of the triangle
                    var cotanVoronoiAlpha = yToYNext.Cotan(-xiToY);
                    var cotanVoronoiBeta = (-yToYNext).Cotan(-xiToYnext);
                    
                    var voronoiArea = 0d; 

                    // Check if triangle is obtuse at X
                    if (Vector3.Dot(xiToY, xiToYnext) < 0d)
                    {
                        // Area of triange / 2
                        voronoiArea += Vector3.Cross(xiToY, xiToYnext).Length() * 0.25d;
                    } else if (Vector3.Dot(yToYNext, -xiToY) < 0d || Vector3.Dot(-yToYNext, -xiToYnext) < 0d)
                    {
                        // Area of triange / 4
                        voronoiArea += Vector3.Cross(xiToY, xiToYnext).Length() * 0.125d;
                    }
                    else
                    {
                        // voronoi(V) = 1/8 * (|u|²*cot(beta) + |v|²*cotan(alpha)
                        voronoiArea += 0.125d * (xiToY.LengthSquared() * cotanVoronoiBeta + xiToYnext.LengthSquared() * cotanVoronoiAlpha);
                    }

                    var cotanMeanCurvAlpha = (-xiToYprev).Cotan(yPrevToY);
                    var cotanMeanCurvBeta = (-xiToYnext).Cotan(-yToYNext);

                    // Inner term of sum for K(x) : 0.5 * (cotan(alpha) + cotan(beta)) * (x - xi)
                    var meanCurvatureScalar = 0.5d * (cotanMeanCurvAlpha + cotanMeanCurvBeta);

                    halfedge.Traits.Cotan = meanCurvatureScalar;
                    halfedge.Traits.VoronoiRegionArea = voronoiArea;
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public static Type SelectedLaplacian { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public static bool LaplacianNormalize { get; set; }


        /// <summary>
        /// 
        /// </summary>
        public static void SerializeBinary<T>(string fileName, T data)
        {
            try
            {
                BinaryFormatter formatter = new BinaryFormatter();
                using (Stream stream = new FileStream(fileName, FileMode.Create, FileAccess.Write, FileShare.None))
                {
                    formatter.Serialize(stream, data);
                    stream.Close();
                }
            }
            catch (Exception ex)
            {
                throw new Exception("DeSerialization Error: " + ex.Message);
            }
        }

        /// <summary>
        /// 
        /// </summary>
        public static void DeSerializeBinary<T>(string fileName, out T data)
        {
            try
            {
                BinaryFormatter formatter = new BinaryFormatter();
                using (Stream stream = new FileStream(fileName, FileMode.Open, FileAccess.Read))
                {
                    data = (T)formatter.Deserialize(stream);
                    stream.Close();
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Serialization Error: " + ex.Message);
            }
        } 
    }   
}
