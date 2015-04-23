

using PacBio.BAM;

public class CheckSWIG
{
   public static void Main()
   {
       var header = new BamHeader();
       header.ToSam();
       System.Console.WriteLine("pbbam SWIG binding to C# worked!");
   }
}
