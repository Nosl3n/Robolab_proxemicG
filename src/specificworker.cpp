/*
 *    Copyright (C) 2024 by YOUR NAME HERE
 *
 *    This file is part of RoboComp
 *
 *    RoboComp is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    RoboComp is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with RoboComp.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "specificworker.h"
#include <cppitertools/filter.hpp>
#include <cppitertools/enumerate.hpp>
#include <cppitertools/sliding_window.hpp>
#include <cppitertools/range.hpp>

/**
* \brief Default constructor
*/
SpecificWorker::SpecificWorker(TuplePrx tprx, bool startup_check) : GenericWorker(tprx)
{
	this->startup_check_flag = startup_check;
}
/**
* \brief Default destructor
*/
SpecificWorker::~SpecificWorker()
{
	std::cout << "Destroying SpecificWorker" << std::endl;
}
bool SpecificWorker::setParams(RoboCompCommonBehavior::ParameterList params)
{
	return true;
}
void SpecificWorker::initialize(int period)  /// se repite una vez
{
	std::cout << "Initialize worker" << std::endl;
	this->Period = period;
	if(this->startup_check_flag)
	{
		this->startup_check();
	}
	else
	{
        // Viewer
        viewer = new AbstractGraphicViewer(this->frame, params.GRID_MAX_DIM);
        viewer->add_robot(params.ROBOT_WIDTH, params.ROBOT_LENGTH, 0, 100, QColor("Blue"));
        viewer->setSceneRect(params.GRID_MAX_DIM);
        viewer->show();

        // Grid
        grid.initialize(params.GRID_MAX_DIM, static_cast<int>(params.TILE_SIZE), &viewer->scene);

        // Lidar thread is created
        read_lidar_th = std::thread(&SpecificWorker::read_lidar,this);
        std::cout << __FUNCTION__ << " Started lidar reader" << std::endl;

        // mouse
        connect(viewer, &AbstractGraphicViewer::new_mouse_coordinates, [this](QPointF p)
        {
            qInfo() << "[MOUSE] New global target arrived:" << p;
//            auto paths = grid.compute_k_paths(Grid::Key{0, 0},
//                                              Grid::Key{p.x(), p.y()},
//                                              1,
//                                              params.MIN_DISTANCE_BETWEEN_PATHS,
//                                              true, false);
//            if(not paths.empty())
//            {
//                draw_path(paths.front(), &viewer->scene);
//                this->lcdNumber_length->display((int) paths.front().size());
//            }
            // get key from point and print key values
            auto &&[success, v] = grid.get_cell(grid.point_to_key(p));
            if(success)
                qInfo() << "Cell key: " << v.id << " free: " << v.free << " visited: " << v.visited << "cost: " << v.cost << "hits: " << v.hits << "misses: " << v.misses;

            //Print key

        });
        connect(viewer, &AbstractGraphicViewer::right_click, [this](QPointF p)
        {
            qInfo() <<  "RIGHT CLICK. Cancelling target";
            draw_path({}, &viewer->scene, true);
            cancel_from_mouse = true;
        });
        if(not params.DISPLAY)
            hide();
		timer.start(params.PERIOD);

        //Reset odometry
        try
        { 
            lidarodometry_proxy->reset();
        }
        catch (const Ice::Exception &e)
        {
            std::cout << "Error reading from LidarOdometry" << e << std::endl;
        }
	}

    grid.id_position_map[1] = std::make_pair(-1000, 1000);
    grid.id_position_map[1] = std::make_pair(1000, 1000);
    grid.id_position_map[1] = std::make_pair(0, 2000);
}
void SpecificWorker::compute() ///solo
{
    /// read LiDAR
    auto res_ = buffer_lidar_data.try_get();
    if (not res_.has_value())  {   /*qWarning() << "No data Lidar";*/ return; }
    auto points = res_.value();
    


    std::pair<Eigen::Transform<double, 3, 1>, Eigen::Transform<double, 3, 1>> robot_pose_and_change;

    if(auto res = get_robot_pose_and_change(); res.has_value())
        robot_pose_and_change = res.value();
    else
    {
        qWarning() << __FUNCTION__ << "No robot pose available. Returning. Check LidarOdometry component status";
        return;
    }

    // /// transform target to robot's frame
    // target = transform_target_to_global_frame(robot_pose_and_change.first, target);    // transform target to robot's frame

    /// clear grid and update it
    mutex_path.lock();
        grid.clear();  // sets all cells to initial values
        grid.update_map(points, Eigen::Vector2f{0.0, 0.0}, params.MAX_LIDAR_RANGE);
        grid.update_costs( params.ROBOT_SEMI_WIDTH, true, robot_pose_and_change);
        grid.contabilizarPosicionActual();
    mutex_path.unlock();

    this->hz = fps.print("FPS:", 3000);
    this->lcdNumber_hz->display(this->hz);
}

///////////////////////////////////////////////////////////////////////////////////////////////
void SpecificWorker::read_lidar()
{
    auto wait_period = std::chrono::milliseconds (this->Period);
    while(true)
    {
        try
        {
            auto data = lidar3d_proxy->getLidarDataWithThreshold2d(params.LIDAR_NAME_LOW,
                                                                   params.MAX_LIDAR_LOW_RANGE,
                                                                   params.LIDAR_LOW_DECIMATION_FACTOR);
            auto data_helios = lidar3d1_proxy->getLidarDataWithThreshold2d(params.LIDAR_NAME_HIGH,
                                                                           params.MAX_LIDAR_HIGH_RANGE,
                                                                           params.LIDAR_HIGH_DECIMATION_FACTOR);
                                                                          

/*             RoboCompLidar3D::TData data;
            RoboCompLidar3D::TData data_helios; */

            // concatenate both lidars
            data.points.insert(data.points.end(), data_helios.points.begin(), data_helios.points.end());
            // compute the period to read the lidar based on the current difference with the lidar period. Use a hysteresis of 2ms
            if (wait_period > std::chrono::milliseconds((long) data.period + 2)) wait_period--;
            else if (wait_period < std::chrono::milliseconds((long) data.period - 2)) wait_period++;
            std::vector<Eigen::Vector3f> eig_data(data.points.size());
            for (const auto &[i, p]: data.points | iter::enumerate)
                eig_data[i] = {p.x, p.y, p.z};
            buffer_lidar_data.put(std::move(eig_data));
        }
        catch (const Ice::Exception &e)
        { std::cout << "Error reading from Lidar3D" << e << std::endl; }
        std::this_thread::sleep_for(wait_period);
    }
} // Thread to read the lidar

//////////////////////////////// Draw ///////////////////////////////////////////////////////
void SpecificWorker::draw_path(const std::vector<Eigen::Vector2f> &path, QGraphicsScene *scene, bool erase_only)
{
    static std::vector<QGraphicsEllipseItem*> points;
    for(auto p : points)
        scene->removeItem(p);
    points.clear();

    if(erase_only) return;

    float s = 100;
    auto color = QColor("green");
    for(const auto &p: path)
    {
        auto ptr = scene->addEllipse(-s/2, -s/2, s, s, QPen(color), QBrush(color));
        ptr->setPos(QPointF(p.x(), p.y()));
        points.push_back(ptr);
    }
}
void SpecificWorker::draw_paths(const std::vector<std::vector<Eigen::Vector2f>> &paths, QGraphicsScene *scene, bool erase_only)
{
    static std::vector<QGraphicsEllipseItem*> points;
    static QColor colors[] = {QColor("cyan"), QColor("blue"), QColor("red"), QColor("orange"), QColor("magenta"), QColor("cyan")};
    for(auto p : points)
        scene->removeItem(p);
    points.clear();

    if(erase_only) return;

    float s = 80;
    for(const auto &[i, path]: paths | iter::enumerate)
    {
        // pick a consecutive color
        auto color = colors[i];
        for(const auto &p: path)
        {
            auto ptr = scene->addEllipse(-s/2.f, -s/2.f, s, s, QPen(color), QBrush(color));
            ptr->setPos(QPointF(p.x(), p.y()));
            points.push_back(ptr);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////
RoboCompGridder::Result SpecificWorker::Gridder_getPaths_unlocked(RoboCompGridder::TPoint source,
                                                         RoboCompGridder::TPoint target,
                                                         int max_paths,
                                                         bool tryClosestFreePoint,
                                                         bool targetIsHuman)
{
    //TODO: improve this method to try to find a path even if the target is not free by using the closest free point
    //TODO: if target is human, set safe area around as free
    RoboCompGridder::Result result;
    std::vector<std::vector<Eigen::Vector2f>> paths;

    auto begin = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    qInfo() << __FUNCTION__ << " New plan request: source [" << source.x << source.y << "], target [" << target.x << target.y << "]"
                            << " max_paths: " << max_paths;

    auto [success, msg, source_key, target_key] =
            grid.validate_source_target(Eigen::Vector2f{source.x, source.y},
                                        source.radius,
                                        Eigen::Vector2f{target.x, target.y},
                                        source.radius);
    if (success)
    {
        //check if is line of sight to target free
        if (grid.is_line_of_sigth_to_target_free(source_key,
                                                 target_key,
                                                 params.ROBOT_SEMI_WIDTH))
        {
            paths.emplace_back(grid.compute_path_line_of_sight(source_key, target_key, params.ROBOT_SEMI_LENGTH));
            if(paths.empty())
                msg = "VLOS path not found";
            else
                msg = "VLOS path";
        }
        else
        {
            paths = grid.compute_k_paths(source_key, target_key,
                                         std::clamp(max_paths, 1, params.NUM_PATHS_TO_SEARCH),
                                         params.MIN_DISTANCE_BETWEEN_PATHS,
                                         tryClosestFreePoint,
                                         targetIsHuman);
            if(paths.empty())
                msg = "Djikstra path not found";
            else
                msg = "Djikstra path";
        }
    }
    result.errorMsg= msg;
    result.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // If not success return result with empty paths, error message and timestamp
    if (not success)
        return result;
    else
    {
        // fill Result with data
        result.paths.resize(paths.size());
        for (const auto &[i, path]: paths | iter::enumerate)
        {
            result.paths[i].resize(path.size());
            for (const auto &[j, point]: path | iter::enumerate)
            {
                result.paths[i][j].x = point.x();
                result.paths[i][j].y = point.y();
            }
        }
        qInfo() << __FUNCTION__ << " " << paths.size() << " paths computed in " <<
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now().time_since_epoch()).count() - begin << " ms" << "Status:" << msg.c_str();
        return result;
    }
}
RoboCompGridder::Result SpecificWorker::Gridder_getPaths(RoboCompGridder::TPoint source,
                                                         RoboCompGridder::TPoint target,
                                                         int max_paths,
                                                         bool tryClosestFreePoint,
                                                         bool targetIsHuman)
{
    //TODO: improve this method to try to find a path even if the target is not free by using the closest free point
    //TODO: if target is human, set safe area around as free
    RoboCompGridder::Result result;
    std::vector<std::vector<Eigen::Vector2f>> paths;

    auto begin = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    qInfo() << __FUNCTION__ << " New plan request: source [" << source.x << source.y << "], target [" << target.x << target.y << "]"
            << " max_paths: " << max_paths;
    mutex_path.lock();
    auto [success, msg, source_key, target_key] =
            grid.validate_source_target(Eigen::Vector2f{source.x, source.y},
                                        source.radius,
                                        Eigen::Vector2f{target.x, target.y},
                                        source.radius);
    if (success)
    {
        //check if is line of sight to target free
        if (grid.is_line_of_sigth_to_target_free(source_key,
                                                 target_key,
                                                 params.ROBOT_SEMI_WIDTH))
        {
            paths.emplace_back(grid.compute_path_line_of_sight(source_key, target_key, params.ROBOT_SEMI_LENGTH));
            if(paths.empty())
                msg = "VLOS path not found";
            else
                msg = "VLOS path";
        }
        else
        {
            paths = grid.compute_k_paths(source_key, target_key,
                                         std::clamp(max_paths, 1, params.NUM_PATHS_TO_SEARCH),
                                         params.MIN_DISTANCE_BETWEEN_PATHS,
                                         tryClosestFreePoint,
                                         targetIsHuman);
            if(paths.empty())
                msg = "Djikstra path not found";
            else
                msg = "Djikstra path";
        }
    }
    mutex_path.unlock();
    result.errorMsg = msg;
    result.timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

    // If not success return result with empty paths, error message and timestamp
    if (not success)
        return result;
    else
    {
        // fill Result with data
        result.paths.resize(paths.size());
        for (const auto &[i, path]: paths | iter::enumerate)
        {
            result.paths[i].resize(path.size());
            for (const auto &[j, point]: path | iter::enumerate)
            {
                result.paths[i][j].x = point.x();
                result.paths[i][j].y = point.y();
            }
        }
        qInfo() << __FUNCTION__ << " " << paths.size() << " paths computed in " <<
                std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::system_clock::now().time_since_epoch()).count() - begin << " ms" << "Status:" << msg.c_str();
        return result;
    }
}
bool SpecificWorker::Gridder_LineOfSightToTarget(RoboCompGridder::TPoint source, RoboCompGridder::TPoint target, float robotRadius)
{
    std::lock_guard<std::mutex> lock(mutex_path);
    auto [success, msg, source_key, target_key] =
            grid.validate_source_target(Eigen::Vector2f{source.x, source.y}, source.radius,
                                        Eigen::Vector2f{target.x, target.y}, target.radius);

    if(success)
        return grid.is_line_of_sigth_to_target_free(source_key, target_key, robotRadius);
    else
        return false;
}
RoboCompGridder::TPoint SpecificWorker::Gridder_getClosestFreePoint(RoboCompGridder::TPoint source)
{
    std::lock_guard<std::mutex> lock(mutex_path);
    if(const auto &p = grid.closest_free({source.x, source.y}); p.has_value())
        return {static_cast<float>(p->x()), static_cast<float>(p->y())};
    else
        return {0, 0};  // non valid closest point  TODO: Change return type so failure to find can be signaled
}
RoboCompGridder::TDimensions SpecificWorker::Gridder_getDimensions()
{
    return {static_cast<float>(params.GRID_MAX_DIM.x()),
            static_cast<float>(params.GRID_MAX_DIM.y()),
            static_cast<float>(params.GRID_MAX_DIM.width()),
            static_cast<float>(params.GRID_MAX_DIM.height())};
}
bool SpecificWorker::Gridder_setGridDimensions(RoboCompGridder::TDimensions dimensions)
{
    qInfo() << __FUNCTION__ << " Setting grid dimensions to [" << dimensions.left << dimensions.top << dimensions.width << dimensions.height << "]";
    params.GRID_MAX_DIM = QRectF(dimensions.left, dimensions.top, dimensions.width, dimensions.height);
    //TODO: update grid, clear and reinitialize
    return true;
}
RoboCompGridder::Result SpecificWorker::Gridder_setLocationAndGetPath(RoboCompGridder::TPoint source, RoboCompGridder::TPoint target, bool setFree, RoboCompGridder::TPoint obstacle)
{
    mutex_path.lock();
    auto source_key = grid.point_to_key(Eigen::Vector2f(target.x,target.y));
    auto obstacle_key = grid.point_to_key(Eigen::Vector2f(obstacle.x,obstacle.y));
    auto submap_copy = grid.copy_submap(obstacle_key, obstacle.radius);

    //print submap copy keys and v.free values
//    for (auto &&[k, v]: submap_copy)
//    {
//        std::cout << "key: " << k.first << " " << k.second << " free: " << v.free << std::endl;
//    }

    //set submap to bool setFree
    grid.set_submap(obstacle_key, obstacle.radius, setFree);
    //get paths
    auto result = Gridder_getPaths_unlocked(source, target, 1, true, true);
    //restore submap
    grid.paste_submap(submap_copy);
    mutex_path.unlock();
    return result;
}
bool SpecificWorker::Gridder_IsPathBlocked(RoboCompGridder::TPath path)
{
    std::lock_guard<std::mutex> lock(mutex_path);
    std::vector<Eigen::Vector2f> path_;
    for(const auto &p: path)
        path_.emplace_back(p.x, p.y);
    return grid.is_path_blocked(path_);
}

//SUBSCRIPTION to setVisualObjects method from VisualElementsPub interface
void SpecificWorker::VisualElementsPub_setVisualObjects(RoboCompVisualElementsPub::TData data)
{

    /*
    std::cout << "Publisher: " << data.publisher << std::endl;
    std::cout << "Timestamp Image: " << data.timestampimage << std::endl;
    std::cout << "Timestamp Generated: " << data.timestampgenerated << std::endl;
    std::cout << "Period: " << data.period << std::endl;
    std::cout << "Objects:" << std::endl;
    */

    grid.id_position_map.clear();

    for (const auto& obj : data.objects) {
        int objectId = obj.id;
        float xPos = 0.0, yPos = 0.0;

        // Buscar las coordenadas x_pos y y_pos en el diccionario TAttributes
        auto itXPos = obj.attributes.find("x_pos");
        auto itYPos = obj.attributes.find("y_pos");

        if (itXPos != obj.attributes.end() && itYPos != obj.attributes.end()) {
            xPos = std::stof(itXPos->second);
            yPos = std::stof(itYPos->second);

            // Almacenar en el diccionario idPositionMap
            grid.id_position_map[objectId] = std::make_pair(xPos, yPos);
        }
    }

    for (const auto& entry : grid.id_position_map) {
        int objectId = entry.first;
        float xPos = entry.second.first;
        float yPos = entry.second.second;

        std::cout << "Object ID: " << objectId << ", Position: (" << xPos << ", " << yPos << ")" << std::endl;
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////
int SpecificWorker::startup_check()
{
    std::cout << "Startup check" << std::endl;
    QTimer::singleShot(200, qApp, SLOT(quit()));
    return 0;
}

std::optional<std::pair<Eigen::Transform<double, 3, 1>, Eigen::Transform<double, 3, 1>>> SpecificWorker::get_robot_pose_and_change()
{
    Eigen::Transform<double, 3, 1> robot_pose;
    Eigen::Transform<double, 3, 1> robot_change;
    try
    {
        //  auto pose = lidarodometry_proxy->getFullPoseMatrix();
        const auto pose_and_change = lidarodometry_proxy->getPoseAndChange();
        const auto &pose = pose_and_change.pose;
        const auto &change = pose_and_change.change;
        robot_pose.matrix() << pose.m00, pose.m01, pose.m02, pose.m03,
                               pose.m10, pose.m11, pose.m12, pose.m13,
                               pose.m20, pose.m21, pose.m22, pose.m23,
                               pose.m30, pose.m31, pose.m32, pose.m33;
        robot_change.matrix() << change.m00, change.m01, change.m02, change.m03,
                                 change.m10, change.m11, change.m12, change.m13,
                                 change.m20, change.m21, change.m22, change.m23,
                                 change.m30, change.m31, change.m32, change.m33;
        // qInfo() << __FUNCTION__ << "Odometry Update";
    }
    catch (const Ice::Exception &e)
    {
        std::cout << "Error reading from LidarOdometry" << e << std::endl;
        return {};
    }
    return std::make_pair(robot_pose, robot_change);
}