using Meshes;
using Meshes.Generic;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Meshes.Utilities;
using SharpDX;
using CSparse.Double;
using System.Diagnostics;


namespace Meshes.Algorithms
{    
    /// <summary>
    /// This class implements various triangle mesh smoothing methods. 
    /// </summary>
    public class MeshSmoothing
    {

        /// TODO_A1 Task 7:
        /// 7.	Implement edge-preserving smoothing (medium)
        /// 
        /// Add an additional smoothing method by adding it to
        /// enum Method, IList MethodCollection and the method Smooth.

        public enum Method
        {
            Explicit_Euler, Implicit_Euler, LS_Optimization, TriangleQuality
        }

        /// <summary>
        /// 
        /// </summary>
        public static readonly IList<string> MethodCollection = new List<string>
        {
            Method.Explicit_Euler.ToString(), 
            Method.Implicit_Euler.ToString(), 
            Method.LS_Optimization.ToString(), 
            Method.TriangleQuality.ToString(),
        };


        /// <summary>
        /// Singleton Instance
        /// </summary>
        public static readonly MeshSmoothing Instance = new MeshSmoothing();

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        public bool Smooth(TriangleMesh mesh)
        {
            switch (this.SelectedMethod)
            {
                default:                    
                case Method.Explicit_Euler:
                    this.ExplicitEulerSmoothing(mesh, this.Iterations);
                    break;
                case Method.Implicit_Euler:
                    this.ImplicitEulerSmoothing(mesh);
                    break;
                case Method.LS_Optimization:
                    this.GlobalLeastSquaresSmoothing(mesh);
                    break;
                case Method.TriangleQuality:
                    this.TriangleQualityOptimization(mesh);
                    break;

            }

            return true;
        }

        /// <summary>
        /// 
        /// </summary>
        public double Lambda { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public int Iterations { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public Method SelectedMethod { get; set; }

        /// <summary>
        /// 
        /// </summary>
        public void Clear()
        {
        }

        

        /// <summary>
        /// Singleton, private c'tor
        /// </summary>
        private MeshSmoothing()
        {
            this.Iterations = 1;
        }

        /// <summary>
        /// Performs global least-squares mesh smoothing with the extended matrix
        /// [ L  ] [x'] =       [ 0 ]
        /// [ Wp ]        [ Wp ][ x ]        
        /// </summary>        
        private void GlobalLeastSquaresSmoothing(TriangleMesh mesh)
        {
            /// output vertex positions (after solving/smoothing)
            double[] xx, xy, xz;
            xx = xy = xz = null;

            /// TODO_A1 Task 5:
            /// 5.	Implement least-squares optimization smoothing (medium)
            /// 
            /// Hint: this function is provided in the documentation
            /// as an example how to use the solver!

            /// update mesh           
            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);
        }

        /// <summary>
        /// performs the Explicit Euler (Taubin) smoothing
        /// </summary>
        private void ExplicitEulerSmoothing(TriangleMesh mesh, int iterationCount)
        {
            /// output vertex positions (after solving/smoothing)
            double[] xx, xy, xz;

            /// get current world positions            
            MeshLaplacian.GetEuclideanCoordinates(mesh, out xx, out xy, out xz);

            var laplacian = MeshLaplacian.CreateLaplacian(mesh, -Lambda, 1f);

            foreach (var currentIteration in Enumerable.Range(0, iterationCount))
            {
                double[] dx, dy, dz;

                MeshLaplacian.ComputeDifferentialCoordinates(laplacian.Compress(), xx, xy, xz, out dx, out dy, out dz);

                xx = dx;
                xy = dy;
                xz = dz;
            }

            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);
        }

        /// <summary>
        /// Performs implicit Euler smoothng (Curvature Flow)
        /// </summary>        
        private void ImplicitEulerSmoothing(TriangleMesh mesh)
        {
            /// output vertex positions (after solving/smoothing)
            double[] xx, xy, xz;
            
            MeshLaplacian.GetEuclideanCoordinates(mesh, out xx, out xy, out xz);

            var solver = QR.Create(MeshLaplacian.CreateLaplacian(mesh, Lambda, 1d).Compress());
            
            foreach (var currentIteration in Enumerable.Range(0, Iterations))
            {
                solver.Solve(xx);
                solver.Solve(xy);
                solver.Solve(xz);
            }

            /// update mesh           
            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        /// <returns></returns>
        private static float CalculateVolume(TriangleMesh mesh)
        {
            var faces = mesh.Faces;

            float v = 0.0f;

            foreach (var f in faces)
            {
                var h1 = f.Halfedge;
                var p1 = h1.FromVertex;
                var p2 = h1.Previous.FromVertex;
                var p3 = h1.Next.FromVertex;
                v += CalculateTriangleVolume(p1, p2, p3);
            }

            return v;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="p3"></param>
        /// <returns></returns>
        private static float CalculateTriangleVolume(TriangleMesh.Vertex p1, TriangleMesh.Vertex p2, TriangleMesh.Vertex p3)
        {
            return 1.0f / 6.0f * (-p3.Traits.Position.X * p2.Traits.Position.Y * p1.Traits.Position.Z + p2.Traits.Position.X * p3.Traits.Position.Y * p1.Traits.Position.Z + p3.Traits.Position.X * p1.Traits.Position.Y * p2.Traits.Position.Z
                - p1.Traits.Position.X * p3.Traits.Position.Y * p2.Traits.Position.Z + p2.Traits.Position.X * p1.Traits.Position.Y * p1.Traits.Position.Z + p2.Traits.Position.X * p3.Traits.Position.Y * p1.Traits.Position.Z);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        private void TriangleQualityOptimization(TriangleMesh mesh)
        {
            double[] xx, xy, xz;
            double[] dx, dy, dz;

            MeshLaplacian.GetEuclideanCoordinates(mesh, out xx, out xy, out xz);

            var cotangensL = MeshLaplacian.CreateCotanLaplacian(mesh, 1d, 0d, true);
            MeshLaplacian.ComputeDifferentialCoordinates(cotangensL.Compress(), xx, xy, xz, out dx, out dy, out dz);

            xx = dx.Concat(xx.Select(x => x * Lambda)).ToArray();
            xy = dy.Concat(xy.Select(x => x * Lambda)).ToArray();
            xz = dz.Concat(xz.Select(x => x * Lambda)).ToArray();
            
            // Use uniform lambda weights for the position constraint and identity weights for the smoothing constraitn
            var optimizationSolver = QR.Create(MeshLaplacian.CreateExtendedUniformLaplacian(mesh, Lambda).Compress());

            optimizationSolver.Solve(xx);
            optimizationSolver.Solve(xy);
            optimizationSolver.Solve(xz);

            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);       
        }
    }    
}
