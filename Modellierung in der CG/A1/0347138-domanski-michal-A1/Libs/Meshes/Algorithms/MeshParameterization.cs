//#define MATLAB

using SharpDX;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Drawing;
using CSparse.Double;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Threading.Tasks;

#if MATLAB
using MatlabWrap;
#endif


namespace Meshes.Algorithms
{
    public enum SolverLibrary
    {
        CSparseDotNet,
        CXSparseDotNet,
        LSQRDotNet,
    }

    /// <summary>
    /// This class implements various triangle mesh parametrization methods. 
    /// </summary>
    public class MeshParameterization
    {
        public const SolverLibrary solver = SolverLibrary.CSparseDotNet;

        public enum Method
        {
            Barycentric, LSCM, DCP, LinearABF
        }

        /// <summary>
        /// 
        /// </summary>
        public static readonly IList<string> MethodCollection = new List<string>
        {
            Method.Barycentric.ToString(),
            Method.LSCM.ToString(),
            Method.DCP.ToString(),
            Method.LinearABF.ToString(),
        };

        /// <summary>
        /// 
        /// </summary>
        public Method SelectedMethod
        {
            get;
            set;
        }


        public int P1Index { get; set; }
        public int P2Index { get; set; }
        public Vector2 P1UV { get; set; }
        public Vector2 P2UV { get; set; }


        /// <summary>
        /// the one an only instance of a singleton class
        /// </summary>
        public static readonly MeshParameterization Instance = new MeshParameterization();

        /// <summary>
        /// private constuctor of a singleton
        /// </summary>
        private MeshParameterization()
        {
            this.SelectedMethod = Method.Barycentric;
            this.P1Index = 153;
            this.P2Index = 19;
            this.P1UV = new Vector2(0.5f, 0);
            this.P2UV = new Vector2(0.5f, 1);

#if MATLAB
            /// initalizes a MALAB session connected to this C# program
            /// if you want to see it in your MATLAB window,
            /// (1) first, start MABLAB *before* starting this program
            /// (2) call in MATLAB: enableservice('AutomationServer',true); (two times?)
            /// (3) start this program, then you can push matrices to MATLAB, see examples below
            Matlab.InitMATLAB(true);
#endif
        }



        /// <summary>
        /// Performs parameterization with the chosen method on the given mesh
        /// The mesh is assumed to be non-closed with proper boundary
        /// </summary>
        /// <param name="meshIn">input mesh, after solving its texture coordinates in vertex traits will be adjusted</param>
        /// <param name="meshOut">an flattened output mesh with only X,Y coordinates set, Z is set to 0</param>
        public void PerformParameterization(TriangleMesh meshin, out TriangleMesh meshout)
        {
            switch (this.SelectedMethod)
            {
                default:
                case Method.Barycentric:
                    BarycentricMapping(meshin, out meshout);
                    break;
                case Method.LSCM:
                    LSCM(meshin, out meshout);
                    break;
                case Method.DCP:
                    DCP(meshin, out meshout);
                    break;
                case Method.LinearABF:
                    LinearABF(meshin, out meshout);
                    break;                                    
            }
        }

        /// <summary>
        /// Barycentric Parameterization
        /// Covers barycentric methods which need a fully defined boundary
        /// A particular method can be chosen by creating an appropriate Laplacian matrix
        /// See also (Floater 2003)
        /// </summary>
        /// <param name="meshin">input mesh, after solving its texture coordinates in vertex traits will be adjusted</param>
        /// <param name="meshout">an flattened output mesh with only X,Y coordinates set, Z is set to 0</param>
        private void BarycentricMapping(TriangleMesh meshin, out TriangleMesh meshout)
        {
            /// init an mesh that serves for output of the 2d parametrized mesh
            meshout = meshin.Copy();
            //meshOut = meshIn;            

            /// get lenghts
            int n = meshout.Vertices.Count;
            int m = meshout.Vertices.Where(x => x.OnBoundary).Count();

            /// right hand side (RHS)
            var bu = new double[n];
            var bv = new double[n];
            var b0 = new double[n];

            /// update mesh positions
            MeshLaplacian.UpdateMesh(meshout, bu, bv, b0, bu, bv);
            MeshLaplacian.UpdateMesh(meshin, bu, bv);
        }

        /// <summary>
        /// The Least-Squares Conformal Mapping method, see (Levy et al. 2002)
        /// Performs linear mapping with free boundary
        /// </summary>
        /// <param name="meshIn">input mesh, after solving its texture coordinates in vertex traits will be adjusted</param>
        /// <param name="meshOut">an flattened output mesh with only X,Y coordinates set, Z is set to 0</param>
        private void LSCM(TriangleMesh meshin, out TriangleMesh meshout)
        {
            /// copy mesh for output
            meshout = meshin.Copy();

            /// get fixed points
            int fv0 = this.P1Index;
            int fv1 = this.P2Index;

            /// provide uv's for fixed 2 points            
            var b = new double[]
            {
                this.P1UV.X, this.P1UV.Y, // u1,v1
                this.P2UV.X, this.P2UV.Y, // u2,v2
            };

            /// get counts
            int n = meshout.Vertices.Count;
            int m = meshout.Faces.Count;

            /// output uv-coordinates
            var bu = new double[n];
            var bv = new double[n];
            var b0 = new double[n];

            var A1 = new TripletMatrix(2 * m, 2 * n - 4, 6 * 2 * m);
            var A2 = new TripletMatrix(2 * m, 4, 4 * 2 * m);


            /// update mesh positions and uv's
            MeshLaplacian.UpdateMesh(meshout, bu, bv, b0, bu, bv);
            MeshLaplacian.UpdateMesh(meshin, bu, bv);
        }

        /// <summary>
        /// The Direct Conformal Parameterization (DCP) method, see (Desbrun et al. 2002)
        /// </summary>
        /// <param name="meshin"></param>
        /// <param name="meshout"></param>
        private void DCP(TriangleMesh meshin, out TriangleMesh meshout)
        {
            /// copy the mesh
            meshout = meshin.Copy();

            /// counters
            int n = meshout.Vertices.Count;
            int m = meshout.Faces.Count;

            /// output uv-coordinates
            var bu = new double[n];
            var bv = new double[n];
            var b0 = new double[n];


            /// update mesh positions and uv's
            MeshLaplacian.UpdateMesh(meshout, bu, bv, b0, bu, bv);
            MeshLaplacian.UpdateMesh(meshin, bu, bv);
        }

        /// <summary>
        /// Linear Angle Based Parameterization (LinABF), see (Zayer et al. 2007)
        /// </summary>
        /// <param name="meshin"></param>
        /// <param name="meshout"></param>
        private void LinearABF(TriangleMesh meshin, out TriangleMesh meshout)
        {
            /// copy the mesh
            meshout = meshin.Copy();

            /// counters
            int n = meshout.Vertices.Count;
            int m = meshout.Faces.Count;

            /// output uv-coordinates
            var bu = new double[n];
            var bv = new double[n];
            var b0 = new double[n];

            /// update mesh positions and uv's
            MeshLaplacian.UpdateMesh(meshout, bu, bv, b0, bu, bv);
            MeshLaplacian.UpdateMesh(meshin, bu, bv);
        }

        /// <summary>
        /// Create a geometry bitmap image from a parameterization (without cutting)
        /// See (Gu et al. 2002)
        /// </summary>
        /// <param name="meshin"></param>
        /// <returns></returns>
        public Bitmap GenerateGeometryImage(TriangleMesh meshin)
        {
            var gimg = new Bitmap(512, 512);

            gimg.Save(string.Format("./../../../../Data/gimg_{0}.png", meshin.FileName), System.Drawing.Imaging.ImageFormat.Png);            
            return gimg;
        }



        /// <summary>
        /// Compute matrix M_t for each triangle
        /// (see slides, lecture 7-2, #9)
        /// </summary>
        private static double[,] ComputeMatrixM(TriangleMesh.Vertex[] vertices)
        {
#if DEBUG
            Debug.Assert(vertices.Length == 3);
#endif
            double[,] M = null;

            return M;
        }

        /// <summary>
        /// Creates a 2x3 matrix from two 1x3 row vectors
        /// </summary>
        /// <param name="r0">1x3 row vector</param>
        /// <param name="r1">1x3 row vector</param>
        /// <returns>2x3 rectangular matrix</returns>
        private static double[,] MatrixFromRowVectors(float[] r0, float[] r1)
        {
#if DEBUG
            Debug.Assert(r0.Length == 3);
            Debug.Assert(r1.Length == 3);
#endif
            var m = new double[2, 3]
            {
                {r0[0], r0[1], r0[2]},
                {r1[0], r1[1], r1[2]},
            };
            return m;
        }

        /// <summary>
        /// Creates a 2x3 matrix from two 1x3 row vectors
        /// </summary>
        /// <param name="r0">1x3 row vector</param>
        /// <param name="r1">1x3 row vector</param>
        /// <returns>2x3 rectangular matrix</returns>
        private static double[,] MatrixFromRowVectors(double[] r0, double[] r1)
        {
#if DEBUG
            Debug.Assert(r0.Length == 3);
            Debug.Assert(r1.Length == 3);
#endif
            var m = new double[2, 3]
            {
                {r0[0], r0[1], r0[2]},
                {r1[0], r1[1], r1[2]},
            };
            return m;
        }

        /// <summary>
        /// Multiply a 2x3 matrix with a 3x1 vector.
        /// </summary>
        /// <param name="m">2x3 matrix</param>
        /// <param name="v">3x1 vector</param>
        /// <returns>2x1 vector</returns>
        private static double[] MatrixMultiply(double[,] m, float[] v)
        {
#if DEBUG
            Debug.Assert(v.Length == 3);
            Debug.Assert(m.Length == 6);
#endif
            double[] value = new double[2]
            {
                m[0, 0] * v[0] + m[0, 1] * v[1] + m[0, 2] * v[2],
                m[1, 0] * v[0] + m[1, 1] * v[1] + m[1, 2] * v[2],
            };
            return value;
        }

        /// <summary>
        /// Multiply a 2x3 matrix with a 3x1 vector.
        /// </summary>
        /// <param name="m">2x3 matrix</param>
        /// <param name="v">3x1 vector</param>
        /// <returns>2x1 vector</returns>
        private static double[] MatrixMultiply(double[,] m, double[] v)
        {
#if DEBUG
            Debug.Assert(v.Length == 3);
            Debug.Assert(m.Length == 6);
#endif
            double[] value = new double[2]
            {
                m[0, 0] * v[0] + m[0, 1] * v[1] + m[0, 2] * v[2],
                m[1, 0] * v[0] + m[1, 1] * v[1] + m[1, 2] * v[2],
            };
            return value;
        }
    }
}
