/*~~~~~~~~~~~
CIS441 Fall 2014
Project1F
Creator: Kyle Rosacker
11/4/14
~~~~~~~~~~~*/

#include <iostream>
#include <sstream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkDataSetWriter.h>




double ceil441(double f)
{
    return ceil(f-0.00001);
}


double floor441(double f)
{
    return floor(f+0.00001);
}

	
double max(double a, double b)
{
	if (b>a)
		return b;
	else
		return a;
}


double min(double a, double b)
{
	if (b<a)
		return b;
	else
		return a;
}


double max(double a, double b, double c)
{
	if (b>a)
		return max(b, c);
	else
		return max(a, c);
}


double min(double a, double b, double c)
{
	if (b<a)
		return min(b, c);
	else
		return min(a, c);
}


vtkImageData * NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
	img->SetDimensions(width, height, 1);
	img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);
	
	return img;
}


void WriteImage(vtkImageData *img, const char *filename)
{
    std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}


class Screen
{
    public:
        unsigned char   *buffer;
        double *zBuffer;
        int width, height;
        
        void setzBuff()
        {
            zBuffer = new double[width*height];
            
            for (int i = 0; i < width * height; i++)
            {
                zBuffer[i] = -1.0;
            }
        }
};


double dotProduct(double v[], double u[], int length)
{
    double product = 0.0;

    for (int i=0; i < length; i++) {
        product += (v[i] * u[i]);
    }
    return product;
}


std::vector<double> crossProduct(double *m1, double *m2)
{
    std::vector<double> product(3);
    product[0] = (m1[1] * m2[2]) - (m1[2] * m2[1]);
    product[1] = (m1[2] * m2[0]) - (m1[0] * m2[2]);
    product[2] = (m1[0] * m2[1]) - (m1[1] * m2[0]);

    return product;
}


struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = +0.6;
         lightDir[1] = 0;
         lightDir[2] = +0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 5.3;
         alpha = 7.5;
         S = 1.0;
    };

    double lightDir[3];     // The direction of the light source
    double Ka;              // The coefficient for ambient lighting.
    double Kd;              // The coefficient for diffuse lighting.
    double Ks;              // The coefficient for specular lighting.
    double alpha;           // The exponent term for specular lighting.
    double S;		    // S for the specular calculation
};


class Matrix
{
  public:
   	 double          A[4][4];

   	 void            TransformPoint(const double *ptIn, double *ptOut);
   	 static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
   	 void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}


class Camera
{
  public:
        double near, far;
        double angle;
        double position[3];
        double focus[3];
        double up[3];

        double O[3];
        double v1[3];
        double v2[3];
        double v3[3];

   	Matrix CameraTransform;
        Matrix ViewTransform;
        Matrix DeviceTransform;
	
	void cTransform()
        {
            
            for (int i = 0; i < 3; i++) {
                O[i] = position[i];
            }
        
            double OminusFocus[3];
            for (int i = 0; i < 3; i++) {
                OminusFocus[i] = O[i] - focus[i];
            }
        
            
            std::vector<double> vector1 = crossProduct(up, OminusFocus);
            for (int i = 0; i < 3; i++) {
                v1[i] = vector1[i];
            }
            double total = 0.0;
            for (int i = 0; i < 3; i++) {
                double temp = v1[i];
                total += temp * temp;
            }
            double v1Square = total;
            for (int i = 0; i < 3; i++) {
                if (fabs(v1Square) < 0.00001)
                    v1[i] = 0;
                else
                    v1[i] = v1[i] / sqrt(v1Square);
            }
            
            
            std::vector<double> vector2 = crossProduct(OminusFocus, v1);
            for (int i = 0; i < 3; i++) {
                v2[i] = vector2[i];
            }
            total = 0.0;
            for (int i = 0; i < 3; i++) {
                double temp = v2[i];
                total += temp * temp;
            }
            double v2Square = total;
            for (int i = 0; i < 3; i++) {
                if (fabs(v2Square) < 0.00001)
                    v2[i] = 0;
                else
                    v2[i] = v2[i] / sqrt(v2Square);
            }
        
            
            for (int i = 0; i < 3; i++) {
                v3[i] = OminusFocus[i];
            }
            total = 0.0;
            for (int i = 0; i < 3; i++) {
                double temp = v3[i];
                total += temp * temp;
            }
            double v3Square = total;
            for (int i = 0; i < 3; i++) {
                if (fabs(v3Square) < 0.00001) {
                    v3[i] = 0;
                }
                else {
                    v3[i] = v3[i] / sqrt(v3Square);
                }
            }
        
            double t[3];
            t[0] = 0 - O[0];
            t[1] = 0 - O[1];
            t[2] = 0 - O[2];
        
            CameraTransform.A[0][0] = v1[0];
            CameraTransform.A[0][1] = v2[0];
            CameraTransform.A[0][2] = v3[0]; 
            CameraTransform.A[0][3] = 0;
            CameraTransform.A[1][0] = v1[1];
            CameraTransform.A[1][1] = v2[1];
            CameraTransform.A[1][2] = v3[1];
            CameraTransform.A[1][3] = 0;
            CameraTransform.A[2][0] = v1[2];
            CameraTransform.A[2][1] = v2[2];
            CameraTransform.A[2][2] = v3[2];
            CameraTransform.A[2][3] = 0;
            CameraTransform.A[3][0] = dotProduct(v1, t, 3);
            CameraTransform.A[3][1] = dotProduct(v2, t, 3);
            CameraTransform.A[3][2] = dotProduct(v3, t, 3);
            CameraTransform.A[3][3] = 1;
            
        }

	void vTransform()
        {
            ViewTransform.A[0][0] = 1 / tan(angle/2);
            ViewTransform.A[0][1] = 0;
            ViewTransform.A[0][2] = 0;
            ViewTransform.A[0][3] = 0;
            ViewTransform.A[1][0] = 0;
            ViewTransform.A[1][1] = 1 / tan(angle/2);
            ViewTransform.A[1][2] = 0;
            ViewTransform.A[1][3] = 0;
            ViewTransform.A[2][0] = 0;
            ViewTransform.A[2][1] = 0;
            ViewTransform.A[2][2] = (far + near) / (far - near);
            ViewTransform.A[2][3] = -1;
            ViewTransform.A[3][0] = 0;
            ViewTransform.A[3][1] = 0;
            ViewTransform.A[3][2] = (2 * far * near) / (far - near);
            ViewTransform.A[3][3] = 0;
                    
        }
    
        void dTransform(Screen s)
        {
            DeviceTransform.A[0][0] = s.width / 2;
            DeviceTransform.A[0][1] = 0;
            DeviceTransform.A[0][2] = 0;
            DeviceTransform.A[0][3] = 0;
            DeviceTransform.A[1][0] = 0;
            DeviceTransform.A[1][1] = s.width / 2;
            DeviceTransform.A[1][2] = 0;
            DeviceTransform.A[1][3] = 0;
            DeviceTransform.A[2][0] = 0;
            DeviceTransform.A[2][1] = 0;
            DeviceTransform.A[2][2] = 1;
            DeviceTransform.A[2][3] = 0;
            DeviceTransform.A[3][0] = s.width / 2;
            DeviceTransform.A[3][1] = s.width / 2;
            DeviceTransform.A[3][2] = 0;
            DeviceTransform.A[3][3] = 1;
                    
        }
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}


double phongShading(LightingParameters &lp, double *viewDirection, double *normal)
{
    double adsShading = 0, diffuse = 0, specular = 0;

    diffuse = fabs(dotProduct(lp.lightDir, normal, 3));

    double initR = 2 * dotProduct(lp.lightDir, normal, 3);
    double R[] = {initR * normal[0] - lp.lightDir[0], initR * normal[1] - lp.lightDir[1], initR * normal[2] - lp.lightDir[2]};
    specular = max(0, lp.S * pow(dotProduct(R, viewDirection, 3), lp.alpha));

    adsShading = lp.Ka + lp.Kd * diffuse + lp.Ks * (specular);

    return adsShading;
}


class GetTriSlope
{
    public:
        double x1, y1, z1, rgb1[3], normal1[3];
        double x2, y2, z2, rgb2[3], normal2[3];
        double slope, intercept;
        bool isVertical, invalidLine;
        double minY, maxY;
    
        GetTriSlope()
        {
        }
    
        GetTriSlope(double X1, double Y1, double Z1, double RGB1[], double Normal1[], double X2, double Y2, double Z2, double RGB2[], double Normal2[])
        {
            x1 = X1;
            y1 = Y1;
            z1 = Z1;

            x2 = X2;
            y2 = Y2;
            z2 = Z2;

            for (int i = 0; i < 3; i++) 
	    {
                rgb1[i] = RGB1[i];
                normal1[i] = Normal1[i];
                rgb2[i] = RGB2[i];
                normal2[i] = Normal2[i];
            }
        
            minY = min(y1, y2);
            maxY = max(y1, y2);
        
            if (x2 - x1 == 0)
            {
                isVertical = true;
                slope = x1;
            }

            else 
	    {
                isVertical = false;
                slope = (y2 - y1) / (x2 - x1);
            }

            intercept = y1 - slope * x1;
            if (y2 - y1 == 0)
                invalidLine = true;
            else
                invalidLine = false;
        }

    	bool yCheck(int y) {
            if (minY <= y and y <= maxY)
                return true;
            else
                return false;
        }

        double y_to_x(double y) {
            if (isVertical == true)
                return slope;
            else {
                if (slope == 0)
                    return 0;
                else
                    return (y - intercept) / slope;
            }
        }

        double interpolation(double v1, double v2, double A, double B, double v) {
            return A + ((v - v1) / (v2 - v1)) * (B - A);
        }
    
        double y_to_z(double y) {
            return interpolation(y1, y2, z1, z2, y);
        }
    
        double y_to_r(double y) {
            return interpolation(y1, y2, rgb1[0], rgb2[0], y);
        }
    
        double y_to_g(double y) {
            return interpolation(y1, y2, rgb1[1], rgb2[1], y);
        }
    
        double y_to_b(double y) {
            return interpolation(y1, y2, rgb1[2], rgb2[2], y);
        }
    
        double y_to_nx(double y) {
            return interpolation(y1, y2, normal1[0], normal2[0], y);
        }
    
        double y_to_ny(double y) {
            return interpolation(y1, y2, normal1[1], normal2[1], y);
        }
    
        double y_to_nz(double y) {
            return interpolation(y1, y2, normal1[2], normal2[2], y);
        }
    
        
};


class Triangle
{
    public:
        double X[3];
        double Y[3];
        double Z[3];
        double colors[3][3];
        double normals[3][3];
        double viewDir[3];
        Screen screen;
        LightingParameters lp;
	
        void getScanline()
        {
            double minX = max(ceil441(min(X[0], X[1], X[2])),0);
            double maxX = min(floor441(max(X[0], X[1], X[2])),screen.width-1);
    
            double minY = max(ceil441(min(Y[0], Y[1], Y[2])),0);
            double maxY = min(floor441(max(Y[0], Y[1], Y[2])),screen.height-1);
    
                    
            GetTriSlope tri1 = GetTriSlope(X[0], Y[0], Z[0], colors[0], normals[0], X[1], Y[1], Z[1], colors[1], normals[1]);
            GetTriSlope tri2 = GetTriSlope(X[1], Y[1], Z[1], colors[1], normals[1], X[2], Y[2], Z[2], colors[2], normals[2]);
            GetTriSlope tri3 = GetTriSlope(X[2], Y[2], Z[2], colors[2], normals[2], X[0], Y[0], Z[0], colors[0], normals[0]);
            
            int offset, offsetZ;
            double getMinX, getMaxX, tempX;
            GetTriSlope leftEnd, rightEnd;

            for (int y = minY; y <= maxY; y++)
            {
                offset = screen.width * 3 * y;
                offsetZ = screen.width * y;
                getMinX = screen.width;
                getMaxX = 0;

                if (!tri1.invalidLine and tri1.yCheck(y))
                {
                    tempX = tri1.y_to_x(y);
                    getMinX = tri1.y_to_x(y);
                    getMaxX = getMinX;
                    leftEnd = tri1;
                    rightEnd = tri1;
                }

                if (!tri2.invalidLine and tri2.yCheck(y))
                {
                    tempX = tri2.y_to_x(y);

                    if (tempX < getMinX)
                    {
                        getMinX = tempX;
                        leftEnd = tri2;
                    }

                    if (tempX > getMaxX)
                    {
                        getMaxX = tempX;
                        rightEnd = tri2;
                    }

                }

                if (!tri3.invalidLine and tri3.yCheck(y))
                {
                    tempX = tri3.y_to_x(y);

                    if (tempX < getMinX)
                    {
                        getMinX = tempX;
                        leftEnd = tri3;
                    }

                    if (tempX > getMaxX)
                    {
                        getMaxX = tempX;
                        rightEnd = tri3;
                    }

                }
        
                double lZ = leftEnd.y_to_z(y);
                double lRGB[] = {leftEnd.y_to_r(y), leftEnd.y_to_g(y), leftEnd.y_to_b(y)};
                double lNormal[] = {leftEnd.y_to_nx(y), leftEnd.y_to_ny(y), leftEnd.y_to_nz(y)};
        
                double rZ = rightEnd.y_to_z(y);
                double rRGB[] = {rightEnd.y_to_r(y), rightEnd.y_to_g(y), rightEnd.y_to_b(y)};
                double rNormal[] = {rightEnd.y_to_nx(y), rightEnd.y_to_ny(y), rightEnd.y_to_nz(y)};
        
                        
                int trueMinX = max(ceil441(getMinX), minX);
                int trueMaxX = min(floor441(getMaxX), maxX);

                for (int x = trueMinX; x <= trueMaxX; x++)
                {
                    double convert;
                    if (getMinX != getMaxX)
                        convert = (((double) x) - getMinX) / (getMaxX - getMinX);
                    else
                        convert = 1.0;
        
                    double z = lZ + convert * (rZ - lZ);
                    if (screen.zBuffer[offsetZ + x] <= z)
                    {
                        double trueNormal[3];
                        trueNormal[0] = lNormal[0] + convert * (rNormal[0] - lNormal[0]);
                        trueNormal[1] = lNormal[1] + convert * (rNormal[1] - lNormal[1]);
                        trueNormal[2] = lNormal[2] + convert * (rNormal[2] - lNormal[2]);
                
                        double getPhong = phongShading(lp, viewDir, trueNormal);
                
                        double shadeRGB[3];
                        shadeRGB[0] = min(1.0, (lRGB[0] + convert * (rRGB[0] - lRGB[0])) * getPhong);
                        shadeRGB[1] = min(1.0, (lRGB[1] + convert * (rRGB[1] - lRGB[1])) * getPhong);
                        shadeRGB[2] = min(1.0, (lRGB[2] + convert * (rRGB[2] - lRGB[2])) * getPhong);
                
                        
                
                        screen.zBuffer[offsetZ + x] = z;
                        screen.buffer[offset + (x * 3)] = (unsigned char) ceil441(255.0 * shadeRGB[0]);
                        screen.buffer[offset + (x * 3) + 1] = (unsigned char) ceil441(255.0 * shadeRGB[1]);
                        screen.buffer[offset + (x * 3) + 2] = (unsigned char) ceil441(255.0 * shadeRGB[2]);
                    }
                }
            }
    
            
        }
};


std::vector<Triangle> GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();
/*
vtkDataSetWriter *writer = vtkDataSetWriter::New();
writer->SetInput(pd);
writer->SetFileName("hrc.vtk");
writer->Write();
 */

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}


void SaveImage(vtkImageData *image, int count)
{
    std::string iName;
    std::ostringstream convert;
    convert << "image" << count;
    iName = convert.str();
    const char * c = iName.c_str();
    WriteImage(image, c);
}


Screen InitializeScreen(unsigned char *buffer, int imgW, int imgH)
{
    Screen screen;
    screen.buffer = buffer;
    screen.width = imgW;
    screen.height = imgH;
    screen.setzBuff();

    
    int npixels = imgW*imgH;
    for (int i = 0 ; i < npixels*3 ; i++) 
        screen.buffer[i] = 0;

    return screen;
}


Triangle Transformer(Triangle t, Matrix m, bool debug)
{
    

    Triangle finalTri;
    for (int i = 0; i < 3; i++) 
    {
        double pIn[4];
        pIn[0] = t.X[i];
        pIn[1] = t.Y[i];
        pIn[2] = t.Z[i];
        pIn[3] = 1;
        double pOut[4];
        m.TransformPoint(pIn, pOut);
        finalTri.X[i] = pOut[0] / pOut[3];
        finalTri.Y[i] = pOut[1] / pOut[3];
        finalTri.Z[i] = pOut[2] / pOut[3];
        
    }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            finalTri.colors[i][j] = t.colors[i][j];
            finalTri.normals[i][j] = t.normals[i][j];
        }
    }
    
    return finalTri;
}


int main()
{

    std::vector<Triangle> triangles = GetTriangles();
    int imgW = 1000;
    int imgH = 1000;
    vtkImageData *image = NewImage(imgW, imgH);
    unsigned char *buffer = (unsigned char *) image->GetScalarPointer(0,0,0);
    Screen s;
    Camera c;
    
    for (int i = 0; i < 1000; i++)
    {
        s = InitializeScreen(buffer, imgW, imgH);
        c = GetCamera(i, 1000);
        c.cTransform();
        c.vTransform();
        c.dTransform(s);
    
        Matrix cvdTrans = Matrix::ComposeMatrices(Matrix::ComposeMatrices(c.CameraTransform, c.ViewTransform), c.DeviceTransform);

    for (int j = 0; j < triangles.size(); j++)
        {
            Triangle t = triangles[j];
            Triangle finalTri = Transformer(t, cvdTrans, false);
            for (int k = 0; k < 3; k++) {
                finalTri.viewDir[k] = c.v3[k];
            }
            finalTri.screen = s;
            finalTri.getScanline();
        }

        SaveImage(image, i);
    
        delete [] s.zBuffer;
    }
}
