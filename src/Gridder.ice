//******************************************************************
// 
//  Generated by RoboCompDSL
//  
//  File name: Gridder.ice
//  Source: Gridder.idsl
//
//******************************************************************
#ifndef ROBOCOMPGRIDDER_ICE
#define ROBOCOMPGRIDDER_ICE
module RoboCompGridder
{
	struct TPoint
	{
		float x;
		float y;
		float radius;
	};
	sequence <TPoint> TPath;
	sequence <TPath> TPaths;
	struct TDimensions
	{
		float left;
		float top;
		float width;
		float height;
	};
	struct Result
	{
		TPaths paths;
		long timestamp;
		string error_msg;
		bool valid;
	};
	interface Gridder
	{
		bool IsPathBlocked (TPath path);
		bool LineOfSightToTarget (TPoint source, TPoint target, float robot_radius);
		TPoint getClosestFreePoint (TPoint source);
		TDimensions getDimensions ();
		Result getPaths (TPoint source, TPoint target, int max_paths, bool try_closest_free_point, bool target_is_human);
		bool setGridDimensions (TDimensions dimensions);
	};
};

#endif