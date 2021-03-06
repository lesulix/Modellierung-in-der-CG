﻿using Meshes;
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
        public enum Method
        {
            Explicit_Euler, Implicit_Euler, LS_Optimization, TriangleQuality, EdgeSmoothing
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
            Method.EdgeSmoothing.ToString(),
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
                case Method.EdgeSmoothing:
                    this.EdgeSmoothingOptimization(mesh);
                    break;

            }

            return true;
        }

        /// <summary>
        /// 
        /// </summary>
        public double Lambda { get; set; }

        public double Mu { get; set; }

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

            var laplacianDeflate = MeshLaplacian.CreateLaplacian(mesh, -Lambda, 1f).Compress();
            var laplacianInflate = MeshLaplacian.CreateLaplacian(mesh, Mu, 1f).Compress();

            foreach (var currentIteration in Enumerable.Range(0, iterationCount))
            {
                double[] dx, dy, dz;

                MeshLaplacian.ComputeDifferentialCoordinates(laplacianDeflate, xx, xy, xz, out dx, out dy, out dz);

                xx = dx;
                xy = dy;
                xz = dz;

                MeshLaplacian.ComputeDifferentialCoordinates(laplacianInflate, xx, xy, xz, out dx, out dy, out dz);

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
            
            solver.Solve(xx);
            solver.Solve(xy);
            solver.Solve(xz);

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
            // Create first n rows of the Linear equation system, L_cot * x
            MeshLaplacian.ComputeDifferentialCoordinates(cotangensL.Compress(), xx, xy, xz, out dx, out dy, out dz);

            // Attach the soft constraints, resulting in the rhs : W_L * f | W_P * V_d, where W_L is the identity
            // f are the gradients computed using cotan weights and W_P are constant lambda weights
            xx = dx.Concat(xx.Select(x => x * Lambda)).ToArray();
            xy = dy.Concat(xy.Select(x => x * Lambda)).ToArray();
            xz = dz.Concat(xz.Select(x => x * Lambda)).ToArray();
            
            // W_L*L | W_P where W_L is the identity, L is the uniform laplacian and W_P are constant lambda weights
            var optimizationSolver = QR.Create(MeshLaplacian.CreateExtendedUniformLaplacian(mesh, Lambda).Compress());

            optimizationSolver.Solve(xx);
            optimizationSolver.Solve(xy);
            optimizationSolver.Solve(xz);

            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);       
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="mesh"></param>
        private void EdgeSmoothingOptimization(TriangleMesh mesh)
        {
            double[] xx, xy, xz;
            double[] dx, dy, dz;

            MeshLaplacian.GetEuclideanCoordinates(mesh, out xx, out xy, out xz);

            // Compute the rhs of the LES : W_L * f | W_P * V_d, where f is null, W_L is the identity W_P are constant lambda weights
            xx = xx.Select(_ => 0d).Concat(xx.Select(x => x * Lambda)).ToArray();
            xy = xy.Select(_ => 0d).Concat(xy.Select(x => x * Lambda)).ToArray();
            xz = xz.Select(_ => 0d).Concat(xz.Select(x => x * Lambda)).ToArray();

            // W_L*L | W_P where W_L is the identity, L is the cotans laplacian and W_P are constant lambda weights
            var optimizationSolver = QR.Create(MeshLaplacian.CreateExtendedCotanLaplacian(mesh, Lambda).Compress());

            optimizationSolver.Solve(xx);
            optimizationSolver.Solve(xy);
            optimizationSolver.Solve(xz);

            MeshLaplacian.UpdateMesh(mesh, xx, xy, xz);
        }
    }    
}
