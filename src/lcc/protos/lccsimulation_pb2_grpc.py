# Generated by the gRPC Python protocol compiler plugin. DO NOT EDIT!
"""Client and server classes corresponding to protobuf-defined services."""
import grpc

import lccsimulation_pb2 as lccsimulation__pb2


class LCCSimulationStub(object):
    """
    Main


    """

    def __init__(self, channel):
        """Constructor.

        Args:
            channel: A grpc.Channel.
        """
        self.Initialize = channel.unary_unary(
                '/lccsimulation.LCCSimulation/Initialize',
                request_serializer=lccsimulation__pb2.Empty.SerializeToString,
                response_deserializer=lccsimulation__pb2.InitData.FromString,
                )
        self.Run = channel.unary_stream(
                '/lccsimulation.LCCSimulation/Run',
                request_serializer=lccsimulation__pb2.Empty.SerializeToString,
                response_deserializer=lccsimulation__pb2.RunData.FromString,
                )
        self.Pause = channel.unary_unary(
                '/lccsimulation.LCCSimulation/Pause',
                request_serializer=lccsimulation__pb2.Empty.SerializeToString,
                response_deserializer=lccsimulation__pb2.ControlSimulationResponse.FromString,
                )
        self.Stop = channel.unary_unary(
                '/lccsimulation.LCCSimulation/Stop',
                request_serializer=lccsimulation__pb2.Empty.SerializeToString,
                response_deserializer=lccsimulation__pb2.ControlSimulationResponse.FromString,
                )


class LCCSimulationServicer(object):
    """
    Main


    """

    def Initialize(self, request, context):
        """Missing associated documentation comment in .proto file."""
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')

    def Run(self, request, context):
        """Missing associated documentation comment in .proto file."""
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')

    def Pause(self, request, context):
        """Missing associated documentation comment in .proto file."""
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')

    def Stop(self, request, context):
        """Missing associated documentation comment in .proto file."""
        context.set_code(grpc.StatusCode.UNIMPLEMENTED)
        context.set_details('Method not implemented!')
        raise NotImplementedError('Method not implemented!')


def add_LCCSimulationServicer_to_server(servicer, server):
    rpc_method_handlers = {
            'Initialize': grpc.unary_unary_rpc_method_handler(
                    servicer.Initialize,
                    request_deserializer=lccsimulation__pb2.Empty.FromString,
                    response_serializer=lccsimulation__pb2.InitData.SerializeToString,
            ),
            'Run': grpc.unary_stream_rpc_method_handler(
                    servicer.Run,
                    request_deserializer=lccsimulation__pb2.Empty.FromString,
                    response_serializer=lccsimulation__pb2.RunData.SerializeToString,
            ),
            'Pause': grpc.unary_unary_rpc_method_handler(
                    servicer.Pause,
                    request_deserializer=lccsimulation__pb2.Empty.FromString,
                    response_serializer=lccsimulation__pb2.ControlSimulationResponse.SerializeToString,
            ),
            'Stop': grpc.unary_unary_rpc_method_handler(
                    servicer.Stop,
                    request_deserializer=lccsimulation__pb2.Empty.FromString,
                    response_serializer=lccsimulation__pb2.ControlSimulationResponse.SerializeToString,
            ),
    }
    generic_handler = grpc.method_handlers_generic_handler(
            'lccsimulation.LCCSimulation', rpc_method_handlers)
    server.add_generic_rpc_handlers((generic_handler,))


 # This class is part of an EXPERIMENTAL API.
class LCCSimulation(object):
    """
    Main


    """

    @staticmethod
    def Initialize(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/lccsimulation.LCCSimulation/Initialize',
            lccsimulation__pb2.Empty.SerializeToString,
            lccsimulation__pb2.InitData.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)

    @staticmethod
    def Run(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_stream(request, target, '/lccsimulation.LCCSimulation/Run',
            lccsimulation__pb2.Empty.SerializeToString,
            lccsimulation__pb2.RunData.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)

    @staticmethod
    def Pause(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/lccsimulation.LCCSimulation/Pause',
            lccsimulation__pb2.Empty.SerializeToString,
            lccsimulation__pb2.ControlSimulationResponse.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)

    @staticmethod
    def Stop(request,
            target,
            options=(),
            channel_credentials=None,
            call_credentials=None,
            insecure=False,
            compression=None,
            wait_for_ready=None,
            timeout=None,
            metadata=None):
        return grpc.experimental.unary_unary(request, target, '/lccsimulation.LCCSimulation/Stop',
            lccsimulation__pb2.Empty.SerializeToString,
            lccsimulation__pb2.ControlSimulationResponse.FromString,
            options, channel_credentials,
            insecure, call_credentials, compression, wait_for_ready, timeout, metadata)
