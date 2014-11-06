using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using CSparse.Double;
using SharpDX;

namespace Meshes.Utilities
{
    public static class Vector3Extensions
    {
        public static double Cotan(this Vector3 lhs, Vector3 rhs)
        {
            return (double) Vector3.Dot(lhs, rhs) / Vector3.Cross(lhs, rhs).Length();
        }
    }

    public static class TripletMatrixExtensions
    {
        public static int RowCount(this TripletMatrix matrix)
        {
            return matrix.m;
        }

        public static int ColumnCount(this TripletMatrix matrix)
        {
            return matrix.n;
        }

        public static int RowAt(this TripletMatrix matrix, int index)
        {
            return matrix.i[index];
        }

        public static int ColumnAt(this TripletMatrix matrix, int index)
        {
            return matrix.p[index];
        }

        public static double ValueAt(this TripletMatrix matrix, int index)
        {
            return matrix.x[index];
        }

        public static TripletMatrix ConcatenateBelow(this TripletMatrix top, TripletMatrix bottom)
        {
            if (top.RowCount() != bottom.RowCount())
                throw new ArgumentException("Top matrix columns do not correspond to bottom matrix columns");

            var result = new TripletMatrix(top.RowCount() + bottom.RowCount(), top.ColumnCount(), top.nzmax + bottom.nzmax);
            Enumerable.Range(0, top.nzmax).Apply(i => 
                result.Entry(top.RowAt(i), top.ColumnAt(i), top.ValueAt(i)));
            Enumerable.Range(0, bottom.nzmax).Apply(i => 
                result.Entry(bottom.RowAt(i) + top.RowCount(), bottom.ColumnAt(i), bottom.ValueAt(i)));

            return result;
        }

        public static void MultiplyBy(this TripletMatrix matrix, double scalar)
        {
            Enumerable.Range(0, matrix.nzmax).Apply(i => 
                matrix.Entry(matrix.RowAt(i), matrix.ColumnAt(i), matrix.ValueAt(i) * scalar));
        }

        public static TripletMatrix SetIdentity(this TripletMatrix matrix)
        {
            if (matrix.n != matrix.m)
                throw new ArgumentException("Matrix is not regular");
        
            Enumerable.Range(0, matrix.n).Apply(i => matrix.Entry(i, i, 1));
            return matrix;
        }
    }
}
